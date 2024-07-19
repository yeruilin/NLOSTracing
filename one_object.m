%%% Tracking one single object at a constant speed using FMCW LiDAR

clear;load("one_object.mat");

%% Mean Filtering
windowsize=10;
fft1=zeros(size(calfft_array1));
fft2=fft1;

for ii=windowsize+1:size(fft1,2)-windowsize
    temp=zeros(size(fft1,1),1);
    for jj=-windowsize:windowsize
        temp=temp+calfft_array1(:,ii+jj);
    end
    temp=temp/(2*windowsize+1);
    fft1(:,ii)=temp;
end

for ii=windowsize+1:size(fft2,2)-windowsize
    temp=zeros(size(fft2,1),1);
    for jj=-windowsize:windowsize
        temp=temp+calfft_array2(:,ii+jj);
    end
    temp=temp/(2*windowsize+1);
    fft2(:,ii)=temp;
end


%% Find Peaks
delta_t=1e-9; %sampling interval
lambda=1550/1e9;%m
c=3e8; %m/s
lambda_gamma=800e-9;%nm/s
fr = 200.323064e6; % optical comb frequency interval
twant_len=t_use;
T=12.5e-3; %chirp period， second
resolution=c/(twant_len*fr); %distance resolution

scanning_num=size(fft1,1);
first_peak1=zeros(scanning_num,1);
first_peak2=zeros(scanning_num,1);
distance_array1=zeros(scanning_num,1);
distance_array2=zeros(scanning_num,1);
peak_indices_array1=zeros(scanning_num,2);
peak_indices_array2=zeros(scanning_num,2);
rate=zeros(scanning_num,1);

for iii=1:scanning_num
    peak_indices1=find_fft_peak(squeeze(fft1(iii,:)),2,100); % 找到一次峰，和另一个多普勒峰
    peak_indices2=find_fft_peak(squeeze(fft2(iii,:)),2,100);
    first_peak1(iii)=(min(peak_indices1))*c/(twant_len*fr)/2; % 一次峰的距离
    first_peak2(iii)=(min(peak_indices2))*c/(twant_len*fr)/2; % 一次峰的距离
    distance_array1(iii)=(max(peak_indices1)-min(peak_indices1))*c/(twant_len*fr)/2; % 非视域峰的距离
    distance_array2(iii)=(max(peak_indices2)-min(peak_indices2))*c/(twant_len*fr)/2; % 非视域峰的距离
    
    peak_indices_array1(iii,1)=min(peak_indices1); peak_indices_array1(iii,2)=max(peak_indices1);
    peak_indices_array2(iii,1)=min(peak_indices2); peak_indices_array2(iii,2)=max(peak_indices2);
    
    df=1/(time_width*delta_t);
    rate(iii)=df*((max(peak_indices2)-min(peak_indices2))-(max(peak_indices1)-min(peak_indices1)))*lambda/4;
end
dis_array=(distance_array2+distance_array1)/2;

% 未修正结果
figure;
subplot(3,2,1);plot(first_peak1);title("first peak distance(upchirp)");xlabel("number");ylabel("distance(m)");
subplot(3,2,2);plot(first_peak2);title("first peak distance(downchirp)");xlabel("number");ylabel("distance(m)");
subplot(3,2,3);plot(distance_array1);title("nlos peak distance(upchirp)");xlabel("number");ylabel("distance(m)");
subplot(3,2,4);plot(distance_array2);title("nlos peak distance(downchirp)");xlabel("number");ylabel("distance(m)");
subplot(3,2,5);plot(dis_array);title("nlos peak distance");xlabel("number");ylabel("distance(m)");
subplot(3,2,6);plot(rate);title("velocity");xlabel("number");ylabel("rate(m/s)");

%% 1D Kalman Filter
A=[1,4*T;0,1]; %define the state matrix
C=eye(2);
var_dist=4e-6;
var_velo=4e-6;calfft_array1=calfft_array1(:,100:end);calfft_array2=calfft_array2(:,100:end);
Q=[var_dist,0;0,var_velo];

R=Q;
Xest=zeros(2,size(dis_array,1)+4);

for j=1:4
    x0 =[dis_array(j,1);rate(j,1)]; %define the initial conditions

    Xest (:,j)=x0;
    P=0.1*eye(2); % inital covariance matrix

    for i=2:size(dis_array,1)/4+1
        P=A*P*A'+Q; %predicting P
        Xest(:,4*(i-1)+j)=A*Xest(:,4*(i-2)+j); %Predicitng the state
        K=P*C'/(C*P*C'+R); %calculating the Kalman gains
        Xest(:,4*(i-1)+j)=Xest(:,4*(i-2)+j)+K*([dis_array(4*(i-2)+j,:)';rate(4*(i-2)+j,:)']-C*Xest(:,4*(i-2)+j));
        P=(eye(2)-K*C)*P;
    end
end

%%% ground truth
gt=zeros(size(dis_array,1),1);
x=1:4:scanning_num;y=dis_array(x);p = polyfit(x, y, 1);gt(x)=p(1)*x+p(2);
x=2:4:scanning_num;y=dis_array(x);p = polyfit(x, y, 1);gt(x)=p(1)*x+p(2);
x=3:4:scanning_num;y=dis_array(x);p = polyfit(x, y, 1);gt(x)=p(1)*x+p(2);
x=4:4:scanning_num;y=dis_array(x);p = polyfit(x, y, 1);gt(x)=p(1)*x+p(2);

figure;plot(dis_array(1:4:scanning_num));hold on;plot(squeeze(Xest(1,5:4:end)));hold on;plot(squeeze(gt(1:4:end)));legend("origin","revised","ground truth");
figure;plot(dis_array(2:4:scanning_num));hold on;plot(squeeze(Xest(1,6:4:end)));hold on;plot(squeeze(gt(2:4:end)));legend("origin","revised","ground truth");
figure;plot(dis_array(3:4:scanning_num));hold on;plot(squeeze(Xest(1,7:4:end)));hold on;plot(squeeze(gt(3:4:end)));legend("origin","revised","ground truth");
figure;plot(dis_array(4:4:scanning_num));hold on;plot(squeeze(Xest(1,8:4:end)));hold on;plot(squeeze(gt(4:4:end)));legend("origin","revised","ground truth");

dis_array=squeeze(Xest(1,5:end));
rate=squeeze(Xest(2,5:end));

%% 3D Positioning
pos1=[-4,4,0]/100;
pos2=[4,4,0]/100;
pos3=[4,-4,0]/100;
pos4=[-4,-4,0]/100;

pos=[pos1; pos2; pos3;pos4];

[V,I]=min(first_peak1(1:4),[],1);
if I==1
    pos=[pos3;pos4;pos1; pos2];
elseif I==2
    pos=[pos2;pos3;pos4;pos1];
elseif I==4
    pos=[pos4;pos1;pos2;pos3];
end

frame_num=scanning_num-3;

% calculate 3D position
dis3d=cal_dis3d(dis_array,pos,frame_num);
dis3d_gt=cal_dis3d(gt,pos,frame_num);

% calculate 3D speed
velo3d=cal_rate3d(rate,dis3d,pos);

%% 3D Kalman Filter
A=[eye(3),T*eye(3);0*eye(3),eye(3)]; %define the state matrix
C=eye(6);
var_dist=0.007;
var_velo=0.004;
Q=[var_dist*eye(3),0*eye(3);0*eye(3),var_velo*eye(3)];

R=Q;% observation error

x0 =[dis3d(1,1);dis3d(1,2);dis3d(1,3);velo3d(1,1);velo3d(1,2);velo3d(1,3);]; %define the initial conditions
Xest=zeros(6,size(dis3d,1)+1);
Xest (:,1)=x0;
P=eye (6); % initial covariance matrix

for i=2:size(Xest,2)
P=A*P*A'+Q; %predicting P
Xest(:,i)=A*Xest(:,i-1); %Predicitng the state
K=P*C'/(C*P*C'+R); %calculating the Kalman gains
Xest(:,i)=Xest(:,i)+K*([dis3d(i-1,:)';velo3d(i-1,:)']-C*Xest(:,i));
P=(eye(6)-K*C)*P;
end

dis3d=Xest(1:3,:)';
velo3d=Xest(4:6,:)';

% average
num=floor(size(dis3d,1)/4);
for i=1:num
    dis3d(i,:)=(dis3d(4*i-3,:)+dis3d(4*i-2,:)+dis3d(4*i-1,:)+dis3d(4*i,:))/4;
    velo3d(i,:)=(velo3d(4*i-3,:)+velo3d(4*i-2,:)+velo3d(4*i-1,:)+velo3d(4*i,:))/4;
    
    dis3d_gt(i,:)=(dis3d_gt(4*i-3,:)+dis3d_gt(4*i-2,:)+dis3d_gt(4*i-1,:)+dis3d_gt(4*i,:))/4;
end

dis3d=dis3d(3:num,:);velo3d=velo3d(3:num,:);dis3d_gt=dis3d_gt(2:num-1,:);

figure;colors=colormap('cool');
colors=colors(end:-1:1,:);colormap(colors);
colorIndices = round(interp1(linspace(0, 1, size(colors, 1)), 1:size(colors, 1), linspace(0, 1, size(dis3d,1)), 'linear', 'extrap'));

dis3d=dis3d*100; % cm
dis3d_gt=dis3d_gt*100;
scatter3(dis3d(:,1),dis3d(:,3),dis3d(:,2), 50, colors(colorIndices, :), 'filled');colorbar;
hold on;
scatter3(dis3d_gt(:,1),dis3d_gt(:,3),dis3d_gt(:,2),'x');
xlabel('X(cm)');ylabel('Z(cm)');zlabel('Y(cm)');% title('Zoomed 3D Trajectory');
zlim([1,3]);ylim([12,14]);
caxis([0.5,1.5]);

legend(["Prediction","Ground truth"]);

dis3d_error=abs(dis3d-dis3d_gt);
fprintf("average 3D position error:%f cm\n",norm(mean(dis3d_error)));
dis3d=dis3d/100;dis3d_gt=dis3d_gt/100; % m

% exportgraphics(gcf,'zoomed_line.pdf','ContentType','vector','BackgroundColor','none');

%% Draw error picture
figure; 
for i=1:3
    plot(10*dis3d_error(:,i));
    
    hold on;
end
xlabel("The frame number");ylabel("Error(mm)");title("3D Position Error");
legend(["\Delta x","\Delta y","\Delta z"]);

real_velo=mean(velo3d,1);

figure; 
for i=1:3
    plot(1e3*(velo3d(:,i)-real_velo(i)));
    
    hold on;
end
xlabel("The frame number");ylabel("Error(mm/s)");title("3D Velocity Error");
legend(["\Delta v_x","\Delta v_y","\Delta v_z"]);
ylim([-1.2,1.2]);