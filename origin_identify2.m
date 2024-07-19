% origin identification based on the last peak
function peak=origin_identify2(peak_indices1,pred_peak1,last_peak)
peak_indices2=peak_indices1;
peaknum=size(peak_indices1,2);
objnum=length(pred_peak1(:));
if peaknum<objnum
    disp(peaknum)
    temp=sort(last_peak);
    dis=temp(2:end)-temp(1:end-1);
    [V,I]=sort(dis);
    for ii=1:objnum-peaknum
        [a,b]=min(abs(temp(I(ii))-peak_indices1),[],2);
        peak_indices2(ii+peaknum)=peak_indices1(b);
    end
    peak_indices2=sort(peak_indices2);
end

[B1,I1]=sort(pred_peak1);
for kk=1:objnum
    peak(I1(kk))=peak_indices2(kk);
end
end