%% createPVscenarioFiles
% Runs the WVM model and puts out the OpenDSS PV scenario files
%
%% Syntax
%  index = createPVscenarioFiles(plantInfoFile,irradianceFile,A_value,circuitFile);
%  index = createPVscenarioFiles();
%
%% Description
% Function to load in the inputs to the WVM (plant info and irradiance
% sensor info), run WVM, create the loadshape file, and the solar scenario
% OpenDSS file.
%
%% Inputs
% * *|plantInfoFile|* - optional input with the link to the MAT file with
% the require PV plant information structure for WVM (see WVM.m and placePVplant.m)
% * *|irradianceFile|* - optional input with the link to the MAT file with
% the require irradiance sensor information structure for WVM (see pvl_WVM.m)
% * *|cloud_speed|* - optional input with a single value of the daily cloud speed 
% * *|circuitFile|* - optional input with the link to the file with the
% OpenDSS circuit
%
%% Outputs
% * *none* - outputs both a .txt loadshape file and a .dss solar scenario
% OpenDSS file
%
%% Copyright 2014
% Georgia Tech Research Corporation, Atlanta, Georgia 30332
% Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
% See the license agreement for full terms and conditions.
%
% Please acknowledge any contributions of the GridPV Toolbox by citing:
% M. J. Reno and K. Coogan, "Grid Integrated Distributed PV (GridPV) Version 2," Sandia National Laboratories SAND2013-20141, 2014.
%
%% Example
% Example run of createPVscenarioFiles
%%
% createPVscenarioFiles('./ExampleCircuit/Ckt24_PV_Central_7_5.mat','./Subfunctions/WVM_subfunctions/Livermore_irr_sensor.mat',10,'.\ExampleCircuit\master_ckt24.dss');
%

function createPVscenarioFiles(plantInfoFile,irradianceFile,cloud_speed,circuitFile)

%% get PV plant info
if ~exist('plantInfoFile')
    [FileName,PathName] = uigetfile('*.mat','Select the file with the PV plant info');
    load([PathName,FileName]);
else
    load(plantInfoFile);
end


%% get sensor data
if ~exist('irradianceFile')
    [FileName,PathName] = uigetfile('*.mat','Select the file with the sensor info');
    load([PathName,FileName]);
else
    load(irradianceFile);
end


%% get A value
if ~exist('cloud_speed')
    cloud_speed = str2double(inputdlg('Insert the cloud speed (m/s)'));
end


%% run WVM
[smooth_irradiance,other_outputs]=pvl_WVM(irr_sensor,plantinfo,cloud_speed);
Power_plant = smooth_irradiance/1000*plantinfo.MW;
Power_plant(isnan(Power_plant)) = 0; % Remove any NaNs caused by WVM's clear-sky model

%% save loadshape
[LoadshapeFileName,LoadshapePathName,FilterIndex] = uiputfile('PVloadshape.txt','Save PV Loadshape file');
if strcmp(plantinfo.powerFactor.type,'fixed')
    csvwrite([LoadshapePathName,LoadshapeFileName],Power_plant);
    loadshapeText = sprintf('new loadshape.PV_Loadshape npts=%i sinterval=%i csvfile="%s" Pbase=%2.2f action=normalize\n',length(Power_plant),round((irr_sensor.time(2)-irr_sensor.time(1))*24*60*60),LoadshapeFileName,plantinfo.MW);
elseif strcmp(plantinfo.powerFactor.type,'VVControl')
    load(plantinfo.powerFactor.filepath);
    csvwrite([LoadshapePathName,LoadshapeFileName],Power_plant);
    loadshapeText = sprintf('new loadshape.PV_Loadshape npts=%i sinterval=%i mult=(file=''%s'') Pbase=%2.2f action=normalize\n',length(Power_plant),round((irr_sensor.time(2)-irr_sensor.time(1))*24*60*60),LoadshapeFileName,plantinfo.MW);
else %power factor as a function of power output or as a function of time
    MVar = makePFprofile(irr_sensor.time,Power_plant,plantinfo.powerFactor.type,plantinfo.powerFactor.filepath,plantinfo.MW);
    csvwrite([LoadshapePathName,LoadshapeFileName(1:end-4),'_P.txt'],Power_plant);
    csvwrite([LoadshapePathName,LoadshapeFileName(1:end-4),'_Q.txt'],-1*MVar);
    loadshapeText = sprintf('new loadshape.PV_Loadshape npts=%i sinterval=%i mult=(file=''%s'') Qmult=(file=''%s'') Pbase=%2.2f Qbase=%2.2f action=normalize\n',length(Power_plant),round((irr_sensor.time(2)-irr_sensor.time(1))*24*60*60),[LoadshapeFileName(1:end-4),'_P.txt'],[LoadshapeFileName(1:end-4),'_Q.txt'],plantinfo.MW,max(MVar));
end

%% get bus locations
location = cd;
[DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;

if ~exist('circuitFile')
    [FileName,PathName] = uigetfile('*.dss','Select the OpenDSS Circuit File of Your Circuit');
    circuitFile = [PathName,FileName];
end
cd(location);
DSSText.command = sprintf('Compile "%s"',circuitFile);
DSSText.command = 'solve';
[busCoordNames busCoordArray] = getBusCoordinatesArray(DSSCircObj);
cd(location);

%% central case
if strcmp(plantinfo.type,'square')
    
    %find closest bus to center of central case
    [junk index] = min(sqrt((busCoordArray(:,2)-plantinfo.polygons(2)).^2+(busCoordArray(:,1)-plantinfo.polygons(1)).^2));
    closestBus = busCoordNames(index);
    Buses = getBusInfo(DSSCircObj,closestBus);
    while Buses.numPhases~=3
        %busCoordArray(:,index) = [];
        busCoordArray(index,:) = [];
        busCoordNames(index) = [];
        [junk index] = min(sqrt((busCoordArray(:,2)-plantinfo.polygons(2)).^2+(busCoordArray(:,1)-plantinfo.polygons(1)).^2));
        closestBus = busCoordNames(index);
        Buses = getBusInfo(DSSCircObj,closestBus);
    end
    
    if strcmp(plantinfo.powerFactor.type,'fixed')
        generatorText = sprintf('new pvsystem.PV bus1=%s irradiance=1 phases=3 kv=%2.2f kVA=%2.2f pf=%2.2f pmpp=%2.2f duty=PV_Loadshape',closestBus{1},Buses.voltage/Buses.voltagePU/1000*sqrt(3),plantinfo.MW*1000*1.1,plantinfo.powerFactor.value,plantinfo.MW*1000); %default to inverter 10% over-rated
    elseif strcmp(plantinfo.powerFactor.type,'VVControl')
        generatorText = sprintf('New XYCurve.myvvc_curve npts=6 Yarray=(%s) XArray=(%s)\n',num2str(VarOutput),num2str(voltage));
        generatorText = [generatorText, sprintf('new pvsystem.PV bus1=%s irradiance=1 phases=3 kv=%2.2f pmpp=%2.2f kva=%2.2f duty=PV_Loadshape\n',closestBus{1},Buses.voltage/Buses.voltagePU/1000*sqrt(3),plantinfo.MW*1000,inverterKVApu*plantinfo.MW*1000)];
        generatorText = [generatorText, sprintf('New InvControl.PVcontrol PVSystemList=pvsystem.PV mode=voltvar VVC_Curve1=myvvc_curve\n')];
        generatorText = [generatorText, 'Set maxcontroliter=500'];
    else %power factor as a function of power output or as a function of time
        generatorText = sprintf('new generator.PV bus1=%s phases=3 kv=%2.2f kw=%2.2f kvar=%2.2f duty=PV_Loadshape',closestBus{1},Buses.voltage/Buses.voltagePU/1000*sqrt(3),plantinfo.MW*1000,max(MVar)*1000);
    end
    
else
%% distributed case
    
    area(:,2) = plantinfo.polygons(:,2);
    area(:,1) = plantinfo.polygons(:,1);
    
    % Get information about each transformer
    Transformers = getTransformerInfo(DSSCircObj);
    
    % Only keep buses inside the area
    IN = inpolygon(busCoordArray(:,1),busCoordArray(:,2),area(:,1),area(:,2));
    busCoordArray = busCoordArray(IN,:);
    busCoordNames = busCoordNames(IN,:);
    
    % Only keep transformers inside the area
    condition = ismember(upper(regexprep({Transformers.bus1},'(\.[0-9]+)','')),upper(busCoordNames));
    Transformers = Transformers(condition);
    
    if length(Transformers)>1 % Allocate PV based on transformer kva rating
        kva = [Transformers.kva];
        totalSystemSize = sum(kva);

        generatorText = [];
        if strcmp(plantinfo.powerFactor.type,'VVControl')
            generatorText = [generatorText, sprintf('New XYCurve.myvvc_curve npts=6 Yarray=(%s) XArray=(%s)\n',num2str(VarOutput),num2str(voltage))];
        end
        for ii=1:length(Transformers)
            Buses = getBusInfo(DSSCircObj,{Transformers(ii).bus2});
            if strcmp(plantinfo.powerFactor.type,'fixed')
                if Transformers(ii).numPhases==1
                    generatorText = [generatorText, sprintf('new pvsystem.PV%s bus1=%s irradiance=1 phases=%1.0f kv=%2.2f kVA=%2.2f pf=%2.2f pmpp=%2.2f duty=PV_Loadshape\n',Transformers(ii).name,Transformers(ii).bus2,Transformers(ii).numPhases,Buses.voltage/Buses.voltagePU/1000,1.1*kva(ii)/totalSystemSize*plantinfo.MW*1000,plantinfo.powerFactor.value,kva(ii)/totalSystemSize*plantinfo.MW*1000)]; %default to inverter 10% over-rated
                else 
                    generatorText = [generatorText, sprintf('new pvsystem.PV%s bus1=%s irradiance=1 phases=%1.0f kv=%2.2f kVA=%2.2f pf=%2.2f pmpp=%2.2f duty=PV_Loadshape\n',Transformers(ii).name,Transformers(ii).bus2,Transformers(ii).numPhases,Buses.voltage/Buses.voltagePU/1000*sqrt(3),1.1*kva(ii)/totalSystemSize*plantinfo.MW*1000,plantinfo.powerFactor.value,kva(ii)/totalSystemSize*plantinfo.MW*1000)]; %default to inverter 10% over-rated
                end
            elseif strcmp(plantinfo.powerFactor.type,'VVControl')
                if Transformers(ii).numPhases==1
                    generatorText = [generatorText, sprintf('new pvsystem.PV%s bus1=%s irradiance=1 phases=%1.0f kv=%2.2f pmpp=%2.2f kva=%2.2f duty=PV_Loadshape\n',Transformers(ii).name,Transformers(ii).bus2,Transformers(ii).numPhases,Buses.voltage/Buses.voltagePU/1000,1.1*kva(ii)/totalSystemSize*plantinfo.MW*1000,kva(ii)/totalSystemSize*plantinfo.MW*1000)];
                else
                    generatorText = [generatorText, sprintf('new pvsystem.PV%s bus1=%s irradiance=1 phases=%1.0f kv=%2.2f pmpp=%2.2f kva=%2.2f duty=PV_Loadshape\n',Transformers(ii).name,Transformers(ii).bus2,Transformers(ii).numPhases,Buses.voltage/Buses.voltagePU/1000*sqrt(3),1.1*kva(ii)/totalSystemSize*plantinfo.MW*1000,kva(ii)/totalSystemSize*plantinfo.MW*1000)];
                end
                generatorText = [generatorText, sprintf('New InvControl.PVcontrol%s PVSystemList=pvsystem.PV%s mode=voltvar VVC_Curve1=myvvc_curve\n',Transformers(ii).name,Transformers(ii).name)];
            else %power factor as a function of power output or as a function of time
                if Transformers(ii).numPhases==1
                    generatorText = [generatorText, sprintf('new generator.PV%s bus1=%s phases=%1.0f kv=%2.2f kw=%2.2f kvar=%2.2f duty=PV_Loadshape\n',Transformers(ii).name,Transformers(ii).bus2,Transformers(ii).numPhases,Buses.voltage/Buses.voltagePU/1000,kva(ii)/totalSystemSize*plantinfo.MW*1000,kva(ii)/totalSystemSize*max(MVar)*1000)];
                else
                    generatorText = [generatorText, sprintf('new generator.PV%s bus1=%s phases=%1.0f kv=%2.2f kw=%2.2f kvar=%2.2f duty=PV_Loadshape\n',Transformers(ii).name,Transformers(ii).bus2,Transformers(ii).numPhases,Buses.voltage/Buses.voltagePU/1000*sqrt(3),kva(ii)/totalSystemSize*plantinfo.MW*1000,kva(ii)/totalSystemSize*max(MVar)*1000)];
                end
            end
        end
    else % Even distribute of PV on load buses
        Loads = getLoadInfo(DSSCircObj);
        
        % Only keep loads inside the area
        condition = ismember(upper(regexprep({Loads.busName},'(\.[0-9]+)','')),upper(busCoordNames));
        Loads = Loads(condition);
        
        generatorText = [];
        if strcmp(plantinfo.powerFactor.type,'VVControl')
            generatorText = [generatorText, sprintf('New XYCurve.myvvc_curve npts=6 Yarray=(%s) XArray=(%s)\n',num2str(VarOutput),num2str(voltage))];
        end
        for ii=1:length(Loads)
            Buses = getBusInfo(DSSCircObj,{Loads(ii).busName});
            if strcmp(plantinfo.powerFactor.type,'fixed')
                if Loads(ii).numPhases==1
                    generatorText = [generatorText, sprintf('new pvsystem.PV%s bus1=%s irradiance=1 phases=%1.0f kv=%2.2f kVA=%2.2f pf=%2.2f pmpp=%2.2f duty=PV_Loadshape\n',Loads(ii).name,Loads(ii).busName,Loads(ii).numPhases,Buses.voltage/Buses.voltagePU/1000,1.1*plantinfo.MW*1000/length(Loads),plantinfo.powerFactor.value,plantinfo.MW*1000/length(Loads))]; %default to inverter 10% over-rated
                else
                    generatorText = [generatorText, sprintf('new pvsystem.PV%s bus1=%s irradiance=1 phases=%1.0f kv=%2.2f kVA=%2.2f pf=%2.2f pmpp=%2.2f duty=PV_Loadshape\n',Loads(ii).name,Loads(ii).busName,Loads(ii).numPhases,Buses.voltage/Buses.voltagePU/1000*sqrt(3),1.1*plantinfo.MW*1000/length(Loads),plantinfo.powerFactor.value,plantinfo.MW*1000/length(Loads))]; %default to inverter 10% over-rated
                end
            elseif strcmp(plantinfo.powerFactor.type,'VVControl')
                if Loads(ii).numPhases==1
                    generatorText = [generatorText, sprintf('new pvsystem.PV%s bus1=%s irradiance=1 phases=%1.0f kv=%2.2f pmpp=%2.2f kva=%2.2f duty=PV_Loadshape\n',Loads(ii).name,Loads(ii).busName,Loads(ii).numPhases,Buses.voltage/Buses.voltagePU/1000,1.1*plantinfo.MW*1000/length(Loads),plantinfo.MW*1000/length(Loads))];
                else
                    generatorText = [generatorText, sprintf('new pvsystem.PV%s bus1=%s irradiance=1 phases=%1.0f kv=%2.2f pmpp=%2.2f kva=%2.2f duty=PV_Loadshape\n',Loads(ii).name,Loads(ii).busName,Loads(ii).numPhases,Buses.voltage/Buses.voltagePU/1000*sqrt(3),1.1*plantinfo.MW*1000/length(Loads),plantinfo.MW*1000/length(Loads))];
                end
                generatorText = [generatorText, sprintf('New InvControl.PVcontrol%s PVSystemList=pvsystem.PV%s mode=voltvar VVC_Curve1=myvvc_curve\n',Loads(ii).name,Loads(ii).name)];
            else %power factor as a function of power output or as a function of time
                if Loads(ii).numPhases==1
                    generatorText = [generatorText, sprintf('new generator.PV%s bus1=%s phases=%1.0f kv=%2.2f kw=%2.2f kvar=%2.2f duty=PV_Loadshape\n',Loads(ii).name,Loads(ii).busName,Loads(ii).numPhases,Buses.voltage/Buses.voltagePU/1000,plantinfo.MW*1000/length(Loads),max(MVar)*1000/length(Loads))];
                else
                    generatorText = [generatorText, sprintf('new generator.PV%s bus1=%s phases=%1.0f kv=%2.2f kw=%2.2f kvar=%2.2f duty=PV_Loadshape\n',Loads(ii).name,Loads(ii).busName,Loads(ii).numPhases,Buses.voltage/Buses.voltagePU/1000*sqrt(3),plantinfo.MW*1000/length(Loads),max(MVar)*1000/length(Loads))];
                end
            end
        end
    end
    
    if strcmp(plantinfo.powerFactor.type,'VVControl')
        generatorText = [generatorText, 'Set maxcontroliter=500'];
    end
    
end

%% create OpenDSS file
[FileName,PathName,FilterIndex] = uiputfile('*.dss','Save the OpenDSS Solar Scenario');
fid = fopen([PathName,FileName],'w');
fprintf(fid,'%s',[loadshapeText,generatorText]);
fclose(fid);

end
