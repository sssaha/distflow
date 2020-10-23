%% PVL_WVM
% Wavelet Variability Model
%
%% Syntax:
% smooth_irradiance = pvl_wvm(irr_sensor,plantinfo,cloud_speed);
%
%% Description
% Computes the spatially-smoothed irradiance to convert a point irradiance sensor to represent a large PV array with decreased ramps.  The method
% uses the wavelet variability model at different time scales to provide all cooresponding smoothing.
%
%% Inputs
% * *|irr_sensor|* is a struct with variables:
% *     irr_sensor.irr: the irradiance measurement
% *     irr_sensor.time: the time stamps (Matlab time vector) for irr_sensor.irr
% *     irr_sensor.Lat: latitude of the sensor
% *     irr_sensor.Lon: longitude of the sensor
% *     irr_sensor.alt: altitude of the sensor
% *     irr_sensor.tilt: tilt angle of the sensor, 0 = flat (e.g., GHI)
% *     irr_sensor.azimuth: azimuth angle of the sensor, 180 = due south
% *     irr.sensor.clear_sky_irradiance:  (optional input) manually enter the clear-sky irradiance (e.g., for an irradiance sensor on a tracking system) 
% *     irr_sensor.UTCoffset=UTC offset
% * *|plantinfo|* is a struct describing the plant to simulate with variables:
% *     plantinfo.tilt: tilt angle of plant modules 
% *     plantinfo.azimuth: azimuth angle of plant modules
% *     plantinfo.clear_sky_irrPOA: (optional input) manually enter the clear-sky irradiance in the module POA (e.g., for tracking systems)
% *     plantinfo.type: 'square','polygon',' or 'discrete'
%       'square' square PV plant with specified number of MWs and PV density
%       'polygon' custom PV plant shape (define vetiticies in lat/lon)
%       'discrete' simulate only certain points (e.g., to replicate output of multiple point sensors)
% *     plantinfo.MW: = MW of PV installed (not necessary for 'discrete' type)
% *     plantinfo.PVdensity: = W PV installed per m2 in plant area (e.g., 41 W/m2 is 1MW per 6 acres) (not necessary for 'discrete' type)
% *     plantinfo.Lat: (only needed for type 'polygon' or 'discrete') latitude of polygon verticies or discrete points
% *     plantinfo.Lon: (only needed for type 'polygon' or 'discrete') longitude of polygon verticies or discrete points
% * *|cloud_speed|* is a single value of the daily cloud speed
%
%% Outputs
% * *|smooth_irradiance*| is the WVM smoothed irradiance representing the average irradiance over the plant footprint. It maintains the time stamps of the input irradiance (irr_sensor.time). 
%
%% BSD License:
% Copyright (c) 2012, The Regents of the University of California
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% Neither the name The Regents of the University of California, the names of its campuses nor any abbreviation thereof, nor the names of the contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%% Example
% This example uses May 18th, 2014 irradiance data collected at Sandia National Laboratories in Livermore, CA to demonstrate use of the wavelet variability model (WVM).
%%
% 
% Livermore=load('.\WVM_subfunctions\Livermore_Sample_GHI.mat');
% irr_sensor.irr=Livermore.GHI; %measured irradiance
% irr_sensor.time=Livermore.dt; %timestamps
% irr_sensor.Lat=37.676208; %sensor latitude
% irr_sensor.Lon=-121.703118; %sensor longitude
% irr_sensor.alt=200; %sensor altitude (in meters)
% irr_sensor.tilt=0; %tilt = 0 for GHI sensor
% irr_sensor.azimuth=180; %180 = due south
% irr_sensor.UTCoffset=-8; %sensor UTC offset
% plantinfo.tilt=37; %assume modules tilted 37 degrees (approximately latitude tilt)
% plantinfo.azimuth=180; %assume modules facing south
% plantinfo.type='square'; %assume a square-shaped PV plant
% plantinfo.MW=30; %assume a 30MW plant
% plantinfo.PVdensity=41; %41 W/m2 = 1MW per 6 acres, which is a standard rule of thumb
% cloud_speed=10; %assume cloud speed of 10 m/s
% [smooth_irradiance,other_outputs]=pvl_WVM(irr_sensor,plantinfo,cloud_speed);
% plot(irr_sensor.time,irr_sensor.irr,'b',irr_sensor.time,smooth_irradiance,'r'); %Zoomed in plot comparing the measured GHI to the WVM output of smoothed POA irradiance.
% legend('measured GHI','WVM smoothed POA');
% set(gca,'xtick',floor(nanmean(irr_sensor.time)):1/(24*12):ceil(nanmean(irr_sensor.time)));
% datetick('x','HH:MM','keepticks','keeplimits');
% xlabel('time of day [HH:MM]');ylabel('Irradiance [W m^{-2}]');
% xlim([floor(nanmean(irr_sensor.time))+10.75/24 floor(nanmean(irr_sensor.time))+11.25/24]);
% title(datestr(nanmean(irr_sensor.time),'mmm-dd-yyyy'));
% 

function [smooth_irradiance,other_outputs]=pvl_WVM(irr_sensor,plantinfo,cloud_speed)

functionPath = which('pvl_WVM');
addpath([functionPath(1:end-9),'WVM_subfunctions\']);
addpath([functionPath(1:end-9),'WVM_subfunctions\nansuite\']);

try
    temp=plantinfo.Lon;
catch
    plantinfo.Lon=irr_sensor.Lon;
    plantinfo.Lat=irr_sensor.Lat;
end

%% compute the clear-sky index
if isfield(irr_sensor,'clear_sky_irradiance')==1
    irr_sensor.tilt=NaN;
    irr_sensor.azimuth=NaN;
    [clear_sky_index]=pvl_WVM_compute_clear_sky_index(irr_sensor.time,irr_sensor.irr,irr_sensor.Lat,irr_sensor.Lon,irr_sensor.alt,irr_sensor.UTCoffset,irr_sensor.tilt,irr_sensor.azimuth,irr_sensor.clear_sky_irradiance);
else
    [clear_sky_index]=pvl_WVM_compute_clear_sky_index(irr_sensor.time,irr_sensor.irr,irr_sensor.Lat,irr_sensor.Lon,irr_sensor.alt,irr_sensor.UTCoffset,irr_sensor.tilt,irr_sensor.azimuth);
end
%% compute the wavelet modes
[wavelet,timeout,tmscales]=pvl_WVM_compute_wavelet(irr_sensor.time,clear_sky_index);

%% compute variablity reduction
try
    dist=otherinputs.dist;
catch
if strcmp(plantinfo.type,'discrete')
    [dist]=pvl_WVM_compute_distances(plantinfo.Lon,plantinfo.Lat,plantinfo.type);
end
if strcmp(plantinfo.type,'square')
    [dist]=pvl_WVM_compute_distances(plantinfo.Lon,plantinfo.Lat,plantinfo.type,plantinfo.MW,plantinfo.PVdensity);
end
if strcmp(plantinfo.type,'polygon')
    [dist]=pvl_WVM_compute_distances(plantinfo.Lon,plantinfo.Lat,plantinfo.type,plantinfo.MW);
end
end

VR=pvl_WVM_compute_VR(dist,tmscales,cloud_speed);

%% smooth wavelets by VR 
for i=1:length(tmscales)
    if i<length(tmscales) %special treatment for last timescale to ensure wavelet modes can be recombined to create simulated wavelet
        wavelet_smooth(i,:)=wavelet(i,:)./sqrt(VR(i));
    else
        wavelet_smooth(i,:)=(wavelet(i,:));
    end
end

%% sum wavelets to creat smoothed clear-sky index
 
[C,ia,ib]=intersect(irr_sensor.time,timeout);

clear_sky_index_smooth=zeros(size(clear_sky_index));
clear_sky_index_smooth(ia)=nansum(wavelet_smooth);

%% compute clear-sky irradiance for plant
irradiance_in=ones(size(irr_sensor.time)).*NaN; %temporary varaible needed in running pvl_WVM_compute_clear_sky_index
if isfield(plantinfo,'clear_sky_irrPOA')==1
    clear_sky_irradiance_smooth=plantinfo.clear_sky_irrPOA;
else
    [clear_sky_indexPOA,clear_sky_irradiance_smooth]=pvl_WVM_compute_clear_sky_index(irr_sensor.time,irradiance_in,mean(plantinfo.Lat),mean(plantinfo.Lon),irr_sensor.alt,irr_sensor.UTCoffset,plantinfo.tilt,plantinfo.azimuth);
end
 
%% combine clear-sky index with clear-sky irradiance to generate smoothed output
try
smooth_irradiance=clear_sky_irradiance_smooth.*clear_sky_index_smooth;
catch
    smooth_irradiance=clear_sky_irradiance_smooth.*clear_sky_index_smooth';
end
%% produce other outputs
other_outputs.clear_sky_index_smooth=clear_sky_index_smooth;
other_outputs.VR=VR;
other_outputs.clear_sky_index=clear_sky_index;
other_outputs.wavelet=wavelet;
other_outputs.wavelet_smooth=wavelet_smooth;
other_outputs.time_temp=timeout;
other_outputs.tmscales=tmscales;
other_outputs.dist=dist;
other_outputs.clear_sky_irr_POA=clear_sky_irradiance_smooth;
other_outputs.time=irr_sensor.time;

