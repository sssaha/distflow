%% DSSStartup
% Function for starting OpenDSS and linking to MATLAB
%
%% Syntax
%  [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
%
%% Description
% Function to start up OpenDSS in the background and bring the program
% handle into MATLAB to allow control of OpenDSS from MATLAB through the
% COM interface.  This function only needs to be executed once per MATLAB
% session.  The same handle to OpenDSS can be used the rest of the session.
% Note: the OpenDSS session started through the COM interface is separate
% from the executable program, so the active circuits and parameters can be
% different between the COM and visual executable.
%
%% Inputs
% * *none*
%
%% Outputs
% * *|DSSCircObj|* is the handle to the object in the OpenDSS program
% containing the circuit object as well as the text object used to the send
% commands to OpenDSS.
% Note: CircuitObj will be empty until Text.command = 'compile example.dss'
% is done to load in an active circuit into the OpenDSS workspace.
% * *|DSSText|* can be used to send commands to OpenDSS through Text.command;
% it can also be called with CircuitObj.Text.command. 
% * *|gridpvPath|* is a string containing the toolbox location
%
%% Example
% Initiating OpenDSS from MATLAB:
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup
%

function [DSSCircObj, DSSText, gridpvPath] = DSSStartup
    
    location = cd;

    %instantiate the DSS Object
    DSSCircObj = actxserver('OpenDSSEngine.DSS');
    DSSCircObj.Allowforms= 0;
    %Start the DSS.   Only needs to be executed the first time w/in a Matlab session
    Start = DSSCircObj.Start(0);
    
    % Check if OpenDSS was started correctly
    if ~Start
        error('OpenDSS startup was unsuccessful.  Please download and install OpenDSS program.')
    end
    
    % Check version number of OpenDSS
    versionNumberArray = cellfun(@str2double,regexp(DSSCircObj.Version,'.^?(\d+)','Tokens'));
    if (versionNumberArray(1)*1e6+versionNumberArray(2)*1e3+versionNumberArray(3)+versionNumberArray(4)/100) < 7006003.31
        error(sprintf('You are running OpenDSS version %i.%i.%i.%i. OpenDSS must be version 7.6.3.31 or higher for some toolbox functionality.  Please install the current version of OpenDSS.',versionNumberArray(1:end-1)))
    end
    
    % Define the text interface
    DSSText = DSSCircObj.Text;
    
    %Clear the circuit from the COM interface to avoid accidentally plotting stale circuits
    DSSText.command = 'clear';
    
    % Get toolbox pathname
    gridpvPath = which('DSSStartup.m');
    gridpvPath = gridpvPath(1:end-12);
    
    cd(location);
end