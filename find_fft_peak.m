function peak_indices=find_fft_peak(data,n,peakdistance)
% data��Ѱ������飬n��Ѱ����Ŀ
peak_indices=zeros(1,n);

new_data=data;
[~,locs]=findpeaks(new_data,'MinPeakDistance',peakdistance); % ����ķ�ǳ�ϸ���������widthЧ�����ã�������distance�ܺ���
peak_values=new_data(locs);
for i=1:n
[val,loc_index]=max(peak_values);
peak_indices(i)=locs(loc_index);
peak_values(loc_index)=-1;
end
peak_indices=sort(peak_indices);% ���к�˳��
end