%% isinterfaceOpenDSS
% Used to check for a valid interface input.
%
%% Syntax
%  isinterface = isinterfaceOpenDSS(DSSCircObj);
%
%% Description
% Used for input parsing. Checks if the input is an OpenDSS COM interface
% and that it is compiled. Returns 1 if it is a compiled OpenDSS object, 0 otherwise.
% If it returns 0, it returns an error indicating whether it failed the
% interface test or the compiled-circuit test.
% 
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
%
%% Outputs
% * *isinterface* - Returns 1 if it is a compiled OpenDSS object, 0 otherwise
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
% Showing interface check
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = 'solve';
% isinterface = isinterfaceOpenDSS(DSSCircObj)
% 

function isinterface = isinterfaceOpenDSS(DSSCircObj)

isinterface = false;

if ispc
    
    % Check if it is a COM interface
    if ~strcmp(class(DSSCircObj),'COM.OpenDSSEngine_DSS')
        error(sprintf('The specified DSSCircObj is not a COM object. \n Please refer to the Manual (DSSStartup -> Outputs) for more information.'))
    else
        % Check for OpenDSS Errors
        if DSSCircObj.Error.Number~=0 && DSSCircObj.AllowForms==0
            warning(sprintf('OpenDSS is reporting an error:\n"%s"',DSSCircObj.Error.Description))
        end
        
        % Check if it is compiled
        if DSSCircObj.NumCircuits == 0
            error(sprintf('The specified OpenDSS circuit object does not contain a compiled circuit. \n You need to run DSSText.command = ''compile YourCircuitName.dss'' \n Please refer to the Manual (DSSStartup -> Outputs) for more information.'))
        else
            % Check if solved
            if DSSCircObj.ActiveCircuit.Solution.Totaliterations == 0
                error(sprintf('The specified OpenDSS circuit object has not been solved.\nTo solve the power flow in OpenDSS, use DSSText.command = ''solve'''))
            else
                % Check if converged solution
                if DSSCircObj.ActiveCircuit.Solution.Converged == 0
                    error(sprintf(['The specified OpenDSS circuit object contains a circuit, but the solution did not converge and is therefore invalid.\n' ...
                    'Try altering the max number of iterations OpenDSS is allowing for solving the power flow.\n Alternatively, there may be an isolated load that is being manually enabled, '...
                    'causing the power flow to not converge.']))
                else
                    isinterface = true;
                    if length(DSSCircObj.ActiveCircuit.Meters.AllNames)==1 && strcmpi(DSSCircObj.ActiveCircuit.Meters.AllNames,'NONE')
                        [msgstr, msgid] = lastwarn;
                        if ~strcmpi(msgid,'GridPV:EnergyMeter') %don't display the warning immediately again
                            warning('GridPV:EnergyMeter','The OpenDSS circuit does not contain an EnergyMeter.  GridPV Toolbox expects an EnergyMeter at the substation.  Refer to OpenDSS documentation for more details.');
                        end
                    end
                end
            end
        end   
    end

end    
