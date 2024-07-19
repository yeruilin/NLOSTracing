%% Find peaks
function peak_indices=find_fft_peak2(data,n,peakdistance,peakheight)
% data:the array to find peaks
% n: peak number
peak_indices=zeros(1,n);

new_data=data;
[~,locs]=findpeaks(new_data,'MinPeakDistance',peakdistance,'MinPeakHeight',peakheight);
peak_values=new_data(locs);
% If the peak number is larger than n, choose the hightest n peaks
if n<=size(locs,2)
    for i=1:n
        [val,loc_index]=max(peak_values);
        peak_indices(i)=locs(loc_index);
        peak_values(loc_index)=-1;
    end
    peak_indices=sort(peak_indices);
else
    peak_indices=sort(squeeze(locs));
end
end