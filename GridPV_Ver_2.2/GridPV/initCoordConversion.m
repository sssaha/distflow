%% initCoordConversion
% Function to initialize the coordinate conversion process
%
%% Syntax
%  initCoordConversion();
%
%% Description
% Function to allow the user to pick between coordinate conversion methods: manual creation or UTM conversion.
%
%% Inputs
% * *None*
%
%% Outputs
% * *None*
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
% 
%%
% initCoordConversion();
%

function initCoordConversion()
    Response = questdlg(sprintf('Is your circuit currently in UTM coordinates?\n\nIf so select ''UTM Conversion.''\n\nOtherwise, select ''Manual Conversion.'''), ...
        'Convert Circuit to GPS', 'UTM Conversion', 'Manual Conversion', 'Manual Conversion');
    switch Response
        case 'UTM Conversion'
            createCircuitCoordConversionUTM();
        case 'Manual Conversion'
            createCircuitCoordConversion();
        otherwise
    end
end