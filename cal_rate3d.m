function velo3d=cal_rate3d(rate,dis3d,pos)
frame_num=size(dis3d,1);
velo3d=zeros(frame_num,3); % 3D speed
cnt=1;

for ii=1:frame_num    
    actual_pos=circshift(pos,[-mod(ii-1,4),0]);
    
    a1=actual_pos(1,:)-dis3d(ii,:);a1=a1/norm(a1);
    a2=actual_pos(2,:)-dis3d(ii,:);a2=a2/norm(a2);
    a3=actual_pos(3,:)-dis3d(ii,:);a3=a3/norm(a3);
    a4=actual_pos(4,:)-dis3d(ii,:);a4=a4/norm(a4);
    a=[a1;a2;a3;a4];
    
    vv=[rate(ii);rate(ii+1);rate(ii+2);rate(ii+3)];
    
    % a'*a has a large condition number, so we use normalization to
    % migrate this problem
    v=inv(a'*a+0.1*eye(3))*(a'*vv);
    
    velo3d(cnt,:)=v;
    cnt=cnt+1;
end
end