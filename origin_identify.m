function peak1=origin_identify(peak_indices1,pred_peak1)
objnum=length(pred_peak1(:));
peak1=zeros(1,objnum);

% If the peak number is equal to object number, match them one by one
if size(peak_indices1,2)==objnum
    [B1,I1]=sort(pred_peak1);
    for kk=1:objnum
        peak1(I1(kk))=peak_indices1(kk);
    end
else
    % otherwise, match the nearest peak
    for kk=1:objnum
        [~,I]=min(abs(pred_peak1(kk)-peak_indices1),[],2);
        peak1(kk)=peak_indices1(I);
    end
end
end