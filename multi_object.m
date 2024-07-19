%%% Tracking multiple objects using FMCW LiDAR

%% Load trajectory
load("multi_object.mat");
objnum=size(velo_gt,1);
frame_num = size(velo_gt,2);

sss=["obj.1","obj.2","obj.3"];

% show ground truth
figure;
for j=1:objnum
    plot3(squeeze(position_gt(j,:,1)), squeeze(position_gt(j,:,3)), squeeze(position_gt(j,:,2)), 'LineWidth', 2);  
    grid on;  
    xlabel('X(cm)');ylabel('Z(cm)');zlabel('Y(cm)');  
    title('Ground Truth'); hold on;
end
legend(sss);

%% Generate FMCW signal
d=20; % cm
pos1=[-d,d,0];
pos2=[d,d,0];
pos3=[d,-d,0];
pos4=[-d,-d,0];

pos=[pos1;pos2;pos3;pos4];

% bin resolution
fmcw_res=0.05; % cm
multi=40;%Doppler multiplier

calfft_array1=rand(frame_num*4,2000);
calfft_array2=rand(frame_num*4,2000);
dis_array_gt=zeros(objnum,frame_num,4);

inten=[1,1,1,1]*20;
% start to simulate
for i=1:frame_num
    for j=1:objnum
        point1=reshape(position_gt(j,i,:),[1,3]);
        rate=squeeze(velo_gt(j,i,:));
        for k=1:4
            index=4*(i-1)+k;
            dis_array_gt(j,i,k)=norm(pos(k,:)-point1);
            ind1=floor(norm(pos(k,:)-point1)/fmcw_res);
            offset1=floor(dot(-rate,point1-pos(k,:))/norm(point1-pos(k,:))*multi);
            calfft_array1(index,ind1-offset1)=calfft_array1(index,ind1-offset1)+inten(j);
            calfft_array2(index,ind1+offset1)=calfft_array2(index,ind1+offset1)+inten(j);
        end
    end
end

%% FMCW NLOS tracking
rate=zeros(objnum,frame_num,4); % radial speed

velo3d=zeros(objnum,frame_num,3);% 3D speed
dis3d=zeros(objnum,frame_num,3); % 3D position

rate(:,1,:)=(distance_array2(:,1,:)-distance_array1(:,1,:))/2/multi;
dis_array=(distance_array1+distance_array2)*fmcw_res/2;

for jj=1:objnum
   dis3d(jj,1,:)=localization(pos,[dis_array(jj,1,1); dis_array(jj,1,2); dis_array(jj,1,3);dis_array(jj,1,4)],[0,0,dis_array(jj,1,1)]);
   velo3d(jj,1,:)=get_velo([rate(jj,1,1);rate(jj,1,2);rate(jj,1,3);rate(jj,1,4)],pos,squeeze(dis3d(jj,1,:)));
end

% Find peak and Origin Identification
lastpeaks1=zeros(4,objnum);lastpeaks2=zeros(4,objnum);
for i=2:frame_num
    % First predict the next state
    pred_peak3d=zeros(objnum,3);
    for j=1:objnum
        pred_peak3d(j,:)=dis3d(j,i-1,:)+deltat*velo3d(j,i-1,:);
    end
    for k=1:4
        % Find peaks
        peak_indices1=find_fft_peak2(squeeze(calfft_array1(4*(i-1)+k,:)),objnum,5,3);
        peak_indices2=find_fft_peak2(squeeze(calfft_array2(4*(i-1)+k,:)),objnum,5,3);
        
        % Predict the peaks at the obtainted signal
        pred_peak1=zeros(1,objnum);pred_peak2=zeros(1,objnum);
        for j=1:objnum
            
            pred_ind=floor(norm(pred_peak3d(j,:)-pos(k,:))/fmcw_res);
            temp=squeeze(pred_peak3d(j,:)-pos(k,:));
            pred_offset=floor(dot(squeeze(velo3d(j,i-1,:)),temp)/norm(temp)*multi);
            pred_peak1(j)=pred_ind-pred_offset;
            pred_peak2(j)=pred_ind+pred_offset;
        end
        % Origin Identification
        if i==2
            peak1=origin_identify(peak_indices1,pred_peak1);lastpeaks1(k,:)=sort(peak1);
            peak2=origin_identify(peak_indices2,pred_peak2);lastpeaks2(k,:)=sort(peak2);
        else
            peak1=origin_identify2(peak_indices1,pred_peak1,squeeze(lastpeaks1(k,:)));lastpeaks1(k,:)=sort(peak1);
            peak2=origin_identify2(peak_indices2,pred_peak2,squeeze(lastpeaks1(k,:)));lastpeaks2(k,:)=sort(peak2);
        end
        for j=1:objnum
            distance_array1(j,i,k)=peak1(j);distance_array2(j,i,k)=peak2(j);
            rate(j,i,k)=(peak2(j)-peak1(j))/2/multi;
            dis_array(j,i,k)=(peak2(j)+peak1(j))*fmcw_res/2;
        end
    end
    % calculate 3D position and speed
    for jj=1:objnum
        xx=localization(pos,[dis_array(jj,i,1); dis_array(jj,i,2); dis_array(jj,i,3);dis_array(jj,i,4)],squeeze(dis3d(jj,i-1,:)));
        dis3d(jj,i,:)=xx;
        
        vv=[rate(jj,i,1);rate(jj,i,2);rate(jj,i,3);rate(jj,i,4)];
        velo3d(jj,i,:)=get_velo(vv,pos,xx);
    end
end

%% Visualization
mycolor=[255,33,222;136,195,102;255,101,0]/256;
colorstr=["spring","summer","autumn","winter"];

figure;
dis_draw=squeeze(dis3d(1,1:20:end,:));colors1=colormap(colorstr(1));colors1=colors1(1:68,:);colorIndices = round(interp1(linspace(0, 1, size(colors1,1)), 1:size(colors1,1), linspace(0, 1, size(dis_draw,1)), 'linear', 'extrap'));p1=scatter3(dis_draw(:,1),dis_draw(:,3),dis_draw(:,2),10,colors1(colorIndices, :),'filled');hold on;
dis_draw=squeeze(dis3d(2,1:20:end,:));colors2=colormap(colorstr(2));colors2=colors2(68:208,:);colorIndices = round(interp1(linspace(0, 1, size(colors2,1)), 1:size(colors2,1), linspace(0, 1, size(dis_draw,1)), 'linear', 'extrap'));p2=scatter3(dis_draw(:,1),dis_draw(:,3),dis_draw(:,2),10,colors2(colorIndices, :),'filled');hold on;
dis_draw=squeeze(dis3d(3,1:20:end,:));colors3=colormap(colorstr(3));colors3=colors3(1:208,:);colorIndices = round(interp1(linspace(0, 1, size(colors3,1)), 1:size(colors3,1), linspace(0, 1, size(dis_draw,1)), 'linear', 'extrap'));p3=scatter3(dis_draw(:,1),dis_draw(:,3),dis_draw(:,2),10,colors3(colorIndices, :),'filled');hold on;

for j=1:objnum
dis_draw=squeeze(dis3d(j,:,:));
plot3(dis_draw(:,1),dis_draw(:,3),dis_draw(:,2),'color',mycolor(j,:));
hold on;
end
xlabel('X(cm)');ylabel('Z(cm)');zlabel('Y(cm)');
legend([p1,p2,p3],sss);

% three point localization:
% pos:position of the scanning points
% dist:radial distance
function xx=localization(pos,dist,lastxx)
    objective = @(x)sum((sqrt(sum((x -pos).^2,2)) - dist).^2);  

    x0=[lastxx(1),lastxx(2),lastxx(3)];
    nonlcon = [];  
    lb = [-inf, -inf, 0];
    ub = [inf, inf, inf];  

    options=optimset('Algorithm','sqp','MaxIter',1e3,'display','off');
    [xx, fval,~,~] = fmincon(objective, x0, [], [], [], [], lb, ub, nonlcon, options);
    disp(fval);
    if fval>0.01
        disp(xx);
    end
end

function velo=get_velo(vv,pos,target_pos)
    target_pos2=reshape(target_pos,size(pos(1,:)));
    a1=pos(1,:)-target_pos2;a1=a1/norm(a1);
    a2=pos(2,:)-target_pos2;a2=a2/norm(a2);
    a3=pos(3,:)-target_pos2;a3=a3/norm(a3);
    a4=pos(4,:)-target_pos2;a4=a4/norm(a4);
    a=[a1;a2;a3;a4];
    velo=-inv(a'*a+0.0*eye(3))*(a'*vv);
end