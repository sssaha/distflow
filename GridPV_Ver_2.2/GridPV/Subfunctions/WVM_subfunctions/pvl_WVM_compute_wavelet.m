function [wavelet,timeout,tmscales]=pvl_WVM_compute_wavelet(timein,clear_sky_index)

%remove NaNs
timein(isnan(clear_sky_index)==1)=[];
clear_sky_index(isnan(clear_sky_index)==1)=[];

%make clear sky index a column vector
sz=size(clear_sky_index);
if sz(1)==1
    clear_sky_index=clear_sky_index';
end

%use symmetric padding to allow for resolution of ends
clear_sky_index2=[flipud(clear_sky_index); clear_sky_index; flipud(clear_sky_index)];

%compute time difference in seconds
tmdiff=round(diff(timein)*60*60*24); 
%compute fraction of values the deviate from the mode
frac_values_diff_from_mode=sum(tmdiff~=mode(tmdiff))./length(tmdiff);

%compute timestep in seconds
secsdiff=mode(tmdiff);

%warn if values deviate from the mode
if frac_values_diff_from_mode>0
warning(['The time sampling interval is not consistent. The mode timestep was ' num2str(secsdiff) ' second(s), but ' num2str(sum(tmdiff~=mode(tmdiff))) ' values (' num2str(round(frac_values_diff_from_mode*100)) '% of values) had different timesteps. The wavelet decomposition may not work as expected.']);
end



%solve 2^j=secsdiff
minj=ceil(log(secsdiff)./log(2));
%only compute wavelets up to 2^12=4096s
lasttimescale=12-minj; 


for j=1:1:lasttimescale
    tmscales(j)=2.^j*secsdiff;
    intervallength=2^j;
    csi_mean(j,:)=moving(clear_sky_index2,intervallength); 
end


wavelet_long=zeros(size(csi_mean)).*NaN;

for i=1:lasttimescale-1
    wavelet_long(i,:)=csi_mean(i,:)-csi_mean(i+1,:);
end
wavelet_long(lasttimescale,:)=csi_mean(lasttimescale,:);
 
for i=1:lasttimescale
    wavelet(i,:)=wavelet_long(i,length(clear_sky_index)+1:2*length(clear_sky_index));
end

timeout=timein;



