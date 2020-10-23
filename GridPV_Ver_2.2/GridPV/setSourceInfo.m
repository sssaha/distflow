function [DSSObj,boolean] = setSourceInfo( DSSObj,sourcename,property,value,varargin )

% This function set the values of properties for a stiff source, not generators.
%
% [DSSObj,boolean] = setSourceInfo( DSSObj,sourcename,property,value,varargin )
% DSSobj= OpenDSS object; 
% sourcename is the name of the sourcess to be changed which can be a cell
% array.
% Property is a string and not case sensitive.
% Value is an array of numbers.
% Length of sourcename and length of value should match, otherwise an error
% is thrown.
% The last input is an optional input, if the user is unsure about the name
% of the source, varargin should be 1 to ensure a sanity check, should be
% avoided if fast operation is required.
% The current implementation does not support multi property setting at the
% same time.
% Current implementation allows to change pu, basekv.
% boolean is set to 1, if property is set properly.
%
% Example:
% setSourceInfo(DSSObj,[{'source'}],'pu',[1.01]) changes
% the pu value of source to 1.01. 

boolean=0;
NameChecker=0;
AdditionalInputLength=length(varargin);
try
    property=lower(property);
catch
    error('String Expected')
end
switch AdditionalInputLength
	case 0
		NameChecker=0;
	case 1
		NameChecker=(varargin{1});
	otherwise
end	
DSSCircuit = DSSObj.ActiveCircuit;
Sources=DSSCircuit.Vsources;
if (NameChecker==1)
	AllSourceNames=Sources.AllNames;
   position=ismember(AllSourceNames,sourcename);
   if (sum(position)~=length(sourcename))
        error('Source Not Found');
    end
end 

if (length(value)~= length(sourcename))
    error ('Data Input Error, number of loads and number of values do not match')
end
for counter= 1:length(value)
    Sources.Name=sourcename{counter};
    
    switch property
        case {'PU','pu'}
            Sources.pu=value(counter);
            boolean=1;
        case 'basekv'
             Sources.BasekV=value(counter);
             boolean=1;
%         case 'kva'
%             Sources.kva=value(counter);
%             boolean=1;
%         case {'pf','PF'}
%             Sources.PF=value(counter);
%             boolean=1;
        otherwise
            warning ('No Property Matched')
    end

end

