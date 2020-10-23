%% IneichenClearSkyModel
% Generates the clear sky irradiance using Ineichen and Perez model 2002
%
%% Syntax
%  GHI = IneichenClearSkyModel(times,latitude,longitude,elevation,Lz);
%
%% Description
% Function to generate the clear sky global horizontal irradiance for a given time period and
% location using the SoDa Linke Turbidity maps
%
%% Inputs
% * *|times|* - matlab datenum (Example: datenum(2011,2,23) ), can be an array of times
% * *|latitude|* - site latitude (decimal degrees)
% * *|longitude|* - site longitude (decimal degrees) (negative for West)
% * *|elevation|* - site elevation (meters)
% * *|Lz|* - standard times zone meridian (120 for PST, 105 for MST, 90 
% CST, and 75 for EST). To find the time zone meridian, just take GMT 
% offset and multiply by -15. (e.g. Eastern time is GMT -5hrs, so the 
% meridian is (-5)*(-15) = 75 degrees.
% * Linke Turbidity images in a folder ('LinkeTurbidity'), images obtained from (http://www.helioclim.org/linke/linke_helioserve.html)
%
%% Outputs
% * *|GHI|* is an array of GHI values for each time in array times
%
%% Example
% Generates the 1-minute GHI profile for Albuquerque for the first week in
% April, 2011.
%%
% times = datenum(2011,4,1):1/(24*60):datenum(2011,4,8);
% GHI=IneichenClearSkyModel(times, 35.04, -106.62, 1617, 105);
% plot(times, GHI,'LineWidth',2); datetick('x','mm/dd','keeplimits','keepticks');
% set(gca,'FontSize',12,'FontWeight','bold');
% ylabel('GHI (W/m^2)','FontSize',12,'FontWeight','bold');
% xlabel('Time','FontSize',12,'FontWeight','bold');
%

function GHI = IneichenClearSkyModel(times,latitude,longitude,elevation,Lz)
%% Parse inputs
p = inputParser; %setup parse structure
p.addRequired('times', @(x) isnumeric(x) && all(x > datenum(1900,1, 1)) && all(x < datenum(2100,1, 1)));
p.addRequired('latitude', @(x) isnumeric(x) && x>-180 && x<180);
p.addRequired('longitude',  @(x) isnumeric(x) && x>-180 && x<180);
p.addRequired('elevation',  @(x) isnumeric(x) && x>=0 && x<9000);
%make sure Lz is one of the 24 time zone meridians
p.addRequired('Lz', @(x) isnumeric(x) && x==165 || x==150 || x==135 || x==120 || x== 105 || x==90 ...
                                      || x==75 || x==60 || x==45 || x==30 || x==15 || x==0 || x==-15 ...
                                      || x==-30 || x==-45 || x==-60 || x==-75 || x==-90 || x==-105 ...
                                      || x==-120 || x==-135 || x==-150 || x==-165 || x==-180);

p.parse(times, latitude, longitude, elevation, Lz); %parse inputs

timeVector = datevec(times);
DOY = ceil( times - datenum(timeVector(:,1)',1,1) );
DOY = DOY';
ToD = timeVector(:,4)+timeVector(:,5)/60+timeVector(:,6)/3600;


%% Extra terrestrial radiation

% From J. W. Spencer, "Fourier series represesnation of the sun," Search, vol. 2, p. 172, 1971.
xx = 2*pi/365*(DOY-1); %in radians
I0=1366.1*(1.00011+0.034221*cos(xx) + 0.00128*sin(xx)-0.000719*cos(2*xx)+0.000077*sin(2*xx)); %solar constant correction for eccentricity and obliquity of Earth's orbit around the sun

% ASCE Standardized Reference
%dr = 1+0.033*cos(2*pi/365 * DOY); % dr = inverse relative distance Earth-Sun (correction for eccentricity of Earth's orbit around the sun) 
%I0=1367.7*dr; %solar constant


%% Zenith Calculation

x = 2*pi/365*(DOY-81); %in radians

phi = pi*latitude/180; %latitude in radians

EoT = 9.87*sin(2*x)-7.53*cos(x)-1.5*sin(x); %minutes
SolarTime = ToD+(Lz+longitude)*4/60+EoT/60;
omega = (SolarTime-12)*15*pi/180;

% Delta From J. W. Spencer, "Fourier series represesnation of the sun," Search, vol. 2, p. 172, 1971.
delta = 0.006918 - 0.399912*cos(xx) + 0.07257*sin(xx) - 0.006758*cos(2*xx) + 0.000907*sin(2*xx) - 0.002697*cos(3*xx) + 0.00148*sin(3*xx);

zenith = acos(sin(phi)*sin(delta) + cos(phi)*cos(delta).*cos(omega) ); %radians

%% Get Linke Turbidity From images for given Lat/Lon
%% Code from Matthew Lave
d=datevec(times);
d1=d(:,2);

Months={'January','February','March','April','May','June','July','August','September','October','November','December'};
La=linspace(90,-90,2160);
Lo=linspace(-180,180,4320);

[C1,I1]=min(abs(latitude-La));
[C2,I2]=min(abs(longitude-Lo));

Tl_months=zeros(12,1);
 
for m=min(d1):1:max(d1)
    pathLocation = which('IneichenClearSkyModel.m');
    pathLocation = pathLocation(1:end-23);
    fname=[pathLocation,'Subfunctions\WVM_subfunctions\LinkeTurbidity\' Months{m} '.tif'];
    Lt=imread(fname);
    L1=Lt(I1,I2);
    L2=double(L1);
    Tl_months(m)=L2./20;
end

% TL=zeros(length(zenith),1);
% for j=1:length(zenith)
%     TL(j)=Tl_months(d1(j));
% end

%interpolate the measured Linke Turbidity
%place the measurement in the middle of each month and interpolate 10 points between each measurement
Tl_months2 = interp1(-0.5:1:14.5, [Tl_months(end-1:end); Tl_months; Tl_months(1:2)], 1:0.1:13,'cubic');
percentageThroughMonth = d(:,3)./eomday(d(:,1),d(:,2));
monthIndex = round((d(:,2)+percentageThroughMonth)*10)-9;
TL = Tl_months2(monthIndex)';

%% Air Mass Calculation
% F. Kasten and A. T. Young, "Revised Optical Air-Mass Tables and Approximation Formula," Applied Optics, vol. 28, pp. 4735-4738, Nov 15 1989.
a=0.50572;
b=6.07995;
c=1.6364;
zenith1=zenith;
zenith1(zenith1>90)=90;
AM=(cos(zenith1)+a.*(90+b-zenith1*180/pi).^(-c)).^(-1);



%% Ineichen and Perez model for GHI
% P. Ineichen and R. Perez, "A new airmass independent formulation for the Linke turbidity coefficient," Solar Energy, vol. 73, pp. 151-157, 2002.
fh1=exp(-elevation/8000);
fh2=exp(-elevation/1250);

cg1=(0.0000509.*elevation+0.868);
cg2=0.0000392.*elevation+0.0387;

GHI=cg1.*I0.*cos(zenith).*exp(-cg2.*AM.*(fh1+fh2*(TL-1))).*exp(0.01.*(AM).^(1.8));
%GHI=cg1.*I0.*cos(zenith).*exp(-cg2.*AM.*(fh1+fh2*(TL*0.65-1))).*exp(0.01.*(AM).^(1.8))*0.95; %changed TL to be 2/3 of original formula for better fit. Also decreased the output just slightly to get closer to a value of one for clearness
GHI(GHI<0)=0;

