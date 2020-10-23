function [clear_sky_index, ClearSkyIrr]=pvl_WVM_compute_clear_sky_index(timein,irradiance_in,latitude,longitude,altitude,UTC,tilt,azimuth,clear_sky_irradiance)

      
 
Time = pvl_maketimestruct(timein, UTC);
Location = pvl_makelocationstruct(latitude, longitude, altitude);

 
[SunAz, SunEl]= pvl_ephemeris(Time, Location);


if exist('clear_sky_irradiance','var')==1
    try
    clear_sky_index=irradiance_in./clear_sky_irradiance;
    catch
            clear_sky_index=irradiance_in./clear_sky_irradiance';
    end
else
    [ClearSkyGHI, ClearSkyDNI, ClearSkyDHI] = pvl_clearsky_ineichen(Time, Location);
    ClearSkyGHI(SunEl<0)=0;
    ClearSkyDNI(SunEl<0)=0;
    ClearSkyDHI(SunEl<0)=0;

    
    if tilt==0
        ClearSkyIrr=ClearSkyGHI;
        try
        clear_sky_index=irradiance_in./ClearSkyGHI;
        catch
            clear_sky_index=irradiance_in./ClearSkyGHI';
        end
    else
        % use Hay/Davies transposition with Ineichen ClearSky Inputs
        doy=pvl_date2doy(Time.year, Time.month, Time.day);
        HExtra = pvl_extraradiation(doy);
        SunZen=90-SunEl;
        SkyDiffuse = pvl_haydavies1980_NaN(tilt,azimuth, ClearSkyDHI, ClearSkyDNI, HExtra, SunZen, SunAz);
        Albedo=0.2;
        GroundRefl = pvl_grounddiffuse(tilt,ClearSkyGHI,Albedo);
        AOI = pvl_getaoi_NaN(tilt,azimuth,SunZen,SunAz);
        Direct = ClearSkyDNI.*cosd(AOI);
        POA=Direct+GroundRefl+SkyDiffuse;
        ClearSkyIrr=POA; 
        try
        clear_sky_index=irradiance_in./POA;
        catch
            clear_sky_index=irradiance_in./POA';
        end
    end
end

clear_sky_index(SunEl<0)=NaN;


%% put limit on clear-sky index when the Sun elevation was < 10 degrees to eliminate very large values
cs1=clear_sky_index(SunEl<10 & SunAz<180);
cs2=clear_sky_index(SunEl<20 & SunEl>10 & SunAz<180);
cs1(cs1>1.2*nanmean(cs2))=1.4*nanmean(cs2);

clear_sky_index(SunEl<10 & SunAz<180)=cs1;

cs3=clear_sky_index(SunEl<10 & SunAz>180);
cs4=clear_sky_index(SunEl<20 & SunEl>10 & SunAz>180);
cs3(cs3>1.2*nanmean(cs4))=1.4*nanmean(cs4);

clear_sky_index(SunEl<10 & SunAz>180)=cs3;
