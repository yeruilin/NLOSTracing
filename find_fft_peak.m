function peak_indices=find_fft_peak(data,n,peakdistance)
% data是寻峰的数组，n是寻峰数目
peak_indices=zeros(1,n);

new_data=data;
[~,locs]=findpeaks(new_data,'MinPeakDistance',peakdistance); % 光梳的峰非常细，因此设置width效果不好，而设置distance很好用
peak_values=new_data(locs);
for i=1:n
[val,loc_index]=max(peak_values);
peak_indices(i)=locs(loc_index);
peak_values(loc_index)=-1;
end
peak_indices=sort(peak_indices);% 排列好顺序
end