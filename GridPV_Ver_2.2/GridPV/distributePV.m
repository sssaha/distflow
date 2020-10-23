%% distributePV
% Allocates PV based off of the load transformer size (kva) 
%
%% Syntax
%  distributePV(totalPVSize,area)
%
%% Description
% Allocates distributed PV spread out around a designated area.  PV is 
% placed at each transformer in the area based off of the load transformer size (kva). 
% The user is asked to select the OpenDSS circuit through the GUI. 
%
%% Inputs
% * *|totalPVSize|* - total size of the distributed PV system in kW
% * *|area|* - matrix of vertices defining the area to distribute the PV inside, 1 row per vertex with [X,Y]
%
%% Outputs
% * text file allocatedPV.txt with the OpenDSS text for PV systems as generators
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
% Distributes the total PV size around the given area.
%%
% area = [1.1732e7 3.708e6; 1.1732e6 3.728e6; 1.1748e7 3.708e6; 1.1748e7 3.728e6];
% totalPVSize = 2e3;
% distributePV(totalPVSize,area);
%

function distributePV(totalPVSize,area)
%% Parse inputs
p = inputParser; %setup parse structure
p.addRequired('totalPVSize', @(x) isnumeric(x) && x>0);
p.addRequired('area', @isnumeric);

p.parse(totalPVSize, area); %parse inputs

%% initiate COM interface (only need to do once when you open MATLAB)
[DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;
    

%% ASsk the user to select the OpenDSS file, then compile it
[FileName,PathName] = uigetfile({'*.txt;*.dss','OpenDSS Files (*.txt,*.dss)'; '*.*',  'All Files (*.*)'},'Select the OpenDSS file with the Circuit');

location = cd;
DSSText.command = sprintf('Compile (%s)',[PathName FileName]);
DSSText.command = 'solve';
cd(location);


%% Get Bus Coordinates
[busCoordNames busCoordArray] = getBusCoordinatesArray(DSSCircObj);


%% Get information about each transformer
Transformers = getTransformerInfo(DSSCircObj);

%% Only keep buses inside the area
IN = inpolygon(busCoordArray(:,1),busCoordArray(:,2),area(:,1),area(:,2));
busCoordArray = busCoordArray(IN,:);
busCoordNames = busCoordNames(IN,:);

%% Only keep transformers inside the area
condition = ismember(regexprep({Transformers.bus1},'(\.[0-9]+)',''),busCoordNames);
Transformers = Transformers(condition);


%% Allocate PV based on transformer kva rating

%kva = mean([[Transformers.kva];repmat(mean([Transformers.kva]),1,length(Transformers))]); %we have some random huge transformers that take all the allocation otherwise
kva = [Transformers.kva];

totalSystemSize = sum(kva);

fid = fopen('allocatedPV.txt','w');
for ii=1:length(Transformers)
    if Transformers(ii).numPhases==1
        fprintf(fid,'new generator.PV%s bus1=%s phases=%1.0f kv=%2.2f kw=%2.2f pf=1 duty=PV_Loadshape\n',Transformers(ii).bus1,Transformers(ii).bus1,Transformers(ii).numPhases,Transformers(ii).bus1Voltage/1000,kva(ii)/totalSystemSize*totalPVSize);
    else
        fprintf(fid,'new generator.PV%s bus1=%s phases=%1.0f kv=%2.2f kw=%2.2f pf=1 duty=PV_Loadshape\n',Transformers(ii).bus1,Transformers(ii).bus1,Transformers(ii).numPhases,Transformers(ii).bus1Voltage/1000*sqrt(3),kva(ii)/totalSystemSize*totalPVSize);
    end
end
fclose(fid);


end