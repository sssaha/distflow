function [TrackingClearSkyIrr,SunEl]=pvl_WVM_compute_single_axis_tracking_clear_sky(timein,latitude,longitude,altitude,UTC,AxisTilt,AxisAzimuth)

     
 
Time = pvl_maketimestruct(timein, UTC);
Location = pvl_makelocationstruct(latitude, longitude, altitude);

 
[ClearSkyGHI, ClearSkyDNI, ClearSkyDHI] = pvl_clearsky_ineichen(Time, Location);
[SunAz, SunEl]= pvl_ephemeris(Time, Location);
ClearSkyGHI(SunEl<0)=0;
ClearSkyDNI(SunEl<0)=0;
ClearSkyDHI(SunEl<0)=0;

        SunZen=90-SunEl;

        MaxAngle=45;
        backtrack=1;
        [TrkrTheta, AOI, SurfTilt, SurfAz,wid] = pvl_singleaxis(SunZen, SunAz, latitude, AxisTilt, AxisAzimuth,MaxAngle,backtrack);



        % use Hay/Davies transposition with Ineichen ClearSky Inputs
        doy=pvl_date2doy(Time.year, Time.month, Time.day);
        HExtra = pvl_extraradiation(doy);
        SkyDiffuse = pvl_haydavies1980_NaN(SurfTilt,SurfAz, ClearSkyDHI, ClearSkyDNI, HExtra, SunZen, SunAz);
        Albedo=0.2;
        GroundRefl = pvl_grounddiffuse(SurfTilt,ClearSkyGHI,Albedo);
        Direct = ClearSkyDNI.*cosd(AOI);
        POA=Direct+GroundRefl+SkyDiffuse;
        TrackingClearSkyIrr=POA; 