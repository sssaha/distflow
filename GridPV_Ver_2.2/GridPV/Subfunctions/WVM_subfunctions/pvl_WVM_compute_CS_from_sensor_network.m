function [cloud_speed]=pvl_WVM_compute_CS_from_sensor_network(timein,Irradiances,Latitudes,Longitudes,Altitudes,UTCoffset,makefigure)
% PVL_WVM_COMPUTE_CS_FROM_SENSOR_NETWORK computes the daily cloud speed
% from one day worth of measurements from an irradiance sensor network.
%
% timein = time vector (use only one timevector for all irradiances)
% Irradiances = matrix whose rows (or columns) are the timeseries of 
%   irradiance measurements
% Latitudes = vector of latitudes of each sensor
% Longitudes = vector of longitudes of each sensor
% Altitudes = vector of altitudes for each sensor
% UTCoffset = UTC offset of timein
% makefigure = (optional input)


sz=size(Irradiances);

if sz(1)>sz(2)
    Irradiances=Irradiances';
end

[sz1 sz2]=size(Irradiances);


%only works for GHI sensors now; will add fix later for POA/tracking
%sensors
tilt=0;
azimuth=180;
tracking=0; 

for i=1:sz1
    try %will fail e.g., if all values at one location are NaNs
    [clear_sky_index_all(i,:)]=pvl_WVM_compute_clear_sky_index(timein,Irradiances(i,:),Latitudes(i),Longitudes(i),Altitudes(i),UTCoffset,tilt,azimuth,tracking);
    [wavelet_all_temp{i},timeout{i},tmscales]=pvl_WVM_compute_wavelet(timein,clear_sky_index_all(i,:));
    %interpolate back to timein (fill in zeros at night)
    [C,ia,ib]=intersect(timeout{i},timeout{1});
    for j=1:length(tmscales)
        if i>1
            wavelet_all{i}(j,:)=zeros(size(wavelet_all{1}(j,:)))*NaN; %fill in NaNs if no data
        end
        wavelet_all{i}(j,ib)=wavelet_all_temp{i}(j,ia);
    end
    end
end

Time = pvl_maketimestruct(timeout{1},UTCoffset);
Location = pvl_makelocationstruct(nanmean(Latitudes),nanmean(Longitudes),nanmean(Altitudes));
[SunAz, SunEl]= pvl_ephemeris(Time, Location);

for i=1:length(tmscales)
    for j=1:sz1
        for k=1:sz1
            try
            cc1=corrcoef(wavelet_all{j}(i,SunEl>10),wavelet_all{k}(i,SunEl>10),'rows','complete');
            corrs{i}(j,k)=cc1(2,1);
            catch
                corrs{i}(j,k)=NaN;
            end
        end
    end
end

[dist]=pvl_WVM_compute_distances(Longitudes,Latitudes,'discrete');

dist2=squareform(dist);
dist2(dist2==0)=NaN;
 
d_t=[];
corrsreshape=[];
for i=1:length(tmscales)
    d_t=[d_t; reshape(dist2,[],1)./tmscales(i)];
    corrsreshape=[corrsreshape; reshape(corrs{i},[],1)];
end

[CS,fval]=fminsearch(@(CS)nanmean((corrsreshape-exp(-d_t./(1/2*CS))).^2),[10]);
cloud_speed=CS;
 
if exist('makefigure','var')
    if makefigure==1
        figure
plot(d_t,corrsreshape,'.','markersize',15);
hold all
plot(((0:0.1:max(d_t))),exp((-(0:0.1:max(d_t)))./(1/2*CS)),'r','linewidth',3);
grid on
xl=get(gca,'xlim');
yl=get(gca,'ylim');
text(0.6*xl(2),-0.2,['CS = ' num2str(round(CS*10)./10) ' m/s'],'fontsize',20,'color','red');
xlabel('distance / timescale [m/s]');
ylabel('correlation')
title([datestr(mean(floor(timein)))]);

 h=axes('position',[0.6 0.55 0.3 0.3]);
 plot(timein,Irradiances)
Time2 = pvl_maketimestruct(timein,UTCoffset);
[SunAz2, SunEl2]= pvl_ephemeris(Time2, Location);
xlim([min(timein(SunEl2>0)) max(timein(SunEl2>0))]);
set(gca,'xtick',floor(mean(timein)):3/24:floor(mean(timein))+1);
datetick('x','HH','keepticks','keeplimits')
ylabel('GHI [W m^{-2}]');
    end
end