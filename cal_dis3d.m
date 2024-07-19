function dis3d=cal_dis3d(dis_array,pos,frame_num)
dis3d=zeros(frame_num,3); % 3D position
cnt=1;
for ii=1:frame_num
    d1=dis_array(ii);
    d2=dis_array(ii+1);
    d3=dis_array(ii+2);
    d4=dis_array(ii+3);
    % initial guess
    x0 = [0, 0, d1];
    
    actual_pos=circshift(pos,[-mod(ii-1,4),0]);

    objective = @(x)sum((sqrt(sum((x -actual_pos).^2,2)) - [d1; d2; d3;d4]).^2);  

    nonlcon = [];  
        
    lb = [-inf, -inf, 0];
    ub = [inf, inf, inf];  
    
    options=optimset('Algorithm','sqp','MaxIter',1e2,'display','off');
    [x, fval,~,~] = fmincon(objective, x0, [], [], [], [], lb, ub, nonlcon, options);
    disp(x);
    disp(fval);
    
    dis3d(cnt,:)=x;
    cnt=cnt+1;
end
end