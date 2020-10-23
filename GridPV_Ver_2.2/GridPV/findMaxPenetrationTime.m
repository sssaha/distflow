%% findMaxPenetrationTime
% Finds the max penetration time
%
%% Syntax
%  index = findMaxPenetrationTime(loadFile,pvFile);
%  index = findMaxPenetrationTime();
%
%% Description
% Function to calculate when the max penetration (PV output / load) time occurs.  User inputs
% the load file and PV output profile, max time is calculated.
%
%% Inputs
% * *|loadFile|* - optional input with the link to the file with the load data
% * *|pvFile|* - optional input with the link to the file with the PV output data
%
%% Outputs
% * *|index|* - the index in the array with the maximum penetration
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
% Finds the maximum penetration time for sample files
%%
% index = findMaxPenetrationTime('ExampleCircuit\LS_ThreePhase.txt','ExampleCircuit\PVloadshape_7_5MW_Central.txt')
%

function index = findMaxPenetrationTime(loadFile,pvFile)

%% Check input arguments and get file paths from user if not in arguments
if nargin==2
    loadData = importdata(loadFile);
    pvOutput = importdata(pvFile);
else
    [FileName,PathName] = uigetfile('*.txt','Select the Loadshape file for the Load');
    loadData = importdata([PathName,FileName]);
    [FileName,PathName] = uigetfile('*.txt','Select the Loadshape file for the Solar Output');
    pvOutput = importdata([PathName,FileName]);
end

%%Get resolutions of loadshape data
pvRes = inputdlg('Enter the resolution of the PV output data in SECONDS:', 'PV Resolution');
pvRes = str2num(pvRes{1})/(24*3600);
loadRes = inputdlg('Enter the resolution of the load data in SECONDS:', 'Load Resolution');
loadRes = str2num(loadRes{1})/(24*3600);

pvTime = pvRes*length(pvOutput);
loadTime = loadRes*length(loadData);

%Interp for the loadshape with lower resolution
if pvRes > loadRes
    pvOutput = interp1(0:pvRes:pvTime-1e-9, pvOutput, 0:loadRes:pvTime-1e-9);
elseif pvRes < loadRes
    loadData = interp1(0:loadRes:loadTime-1e-9, loadData, 0:pvRes:loadTime-1e-9);
end

%Truncate the longer data file to be the same size as the shorter
if length(pvOutput) > length(loadData)
    pvOutput = pvOutput(1:length(loadData));
elseif length(loadData) > length(pvOutput)
    loadData = loadData(1:length(pvOutput));
end

%% Find max time
solarPenetration = reshape(pvOutput,1,[])./reshape(loadData,1,[]);
[maxPenetration index] = max(solarPenetration);
    
end
