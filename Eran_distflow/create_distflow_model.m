function [Bus,Branch,BusVoltage] = create_distflow_model(DSSObj,Baseparam)

%% Retrieving the parameters from OpenDSS Solution, running this code requires to have the GridPV library
BusInfo=getBusInfo(DSSObj);
BusVoltage = [BusInfo.phaseVoltagesPU];
BusVoltage = BusVoltage(BusVoltage>0);
LineInfo=getLineInfo(DSSObj);
BusNames={BusInfo.name};
FromBusName={LineInfo.bus1Name};
ToBusName={LineInfo.bus2Name};
CapConversionFactor=1i*(1e-9)*2*pi*DSSObj.ActiveCircuit.Solution.Frequency;
DSSObj.ActiveCircuit.Vsources.First;
slack_voltage=DSSObj.ActiveCircuit.Vsources.pu;
% LineParamModulator=1000;
sbase=Baseparam.sbase;
zbase=Baseparam.zbase;
delta_index=Baseparam.delta_index;

%% Getting the From and To Node list
f=zeros(1,length(FromBusName));
t=zeros(1,length(FromBusName));

% Giving bus Names to a ID 
BusNametoId = struct();

for k = 1:length(BusInfo)
    disp (BusInfo(k).name)
    BusNametoId.(['b' BusInfo(k).name]) = BusInfo(k).id; 
end

for i = 1:length(f)
    f(i)=BusNametoId.(['b' LineInfo(i).bus1Name]);
    t(i)=BusNametoId.(['b' LineInfo(i).bus2Name]);
end

%% Doing the Node Rebelling
f=f';
t=t';
[nmap,rmap,fnew,tnew]=NodeRelabling(f,t,1,'plots');
[~,I] = sort(tnew);
fnew  = fnew(I);
tnew  = tnew(I);

%%
figure
G=graph(FromBusName,ToBusName);
plot(G,'Layout', 'layered','Sources','sourcebus');
labels = cell(length(BusNames),1);
for k = 1:length(BusNames)
    labels{k} = sprintf('%s<->%d<->%d',G.Nodes.Name{k}, BusNametoId.(['b' G.Nodes.Name{k}]),nmap(BusNametoId.(['b' G.Nodes.Name{k}])));
end
plot(G,'Layout', 'layered', 'Sources', 'sourcebus', 'NodeLabel', labels)
title('Complete Labeling')
set(gca,'xtick',[])
set(gca,'xtick',[], 'xticklabel',[], 'ytick',[], 'yticklabel',[])

%% Getting the To array according to the new mapping 
newf=nmap(f);
newt=nmap(t);

%%
% Each Line is named according to its to node
NewLineNames=strcat('l_',string(newt-1));
% Forming the structure for line R, X , Y and phase information
LineR=struct();
LineX=struct();
LineB=struct();
LineY=struct();
LineP=struct();
LineZ=struct();
LineFlow=struct();
% This is not a very elegant way, but will work for now
for i = 1:length(NewLineNames)
    LineR.(NewLineNames{i})=vec2mat(LineInfo(i).Rmatrix,LineInfo(i).numPhases)* LineInfo(i).length;
    LineX.(NewLineNames{i})=vec2mat(LineInfo(i).Xmatrix,LineInfo(i).numPhases)* LineInfo(i).length;
    LineB.(NewLineNames{i})=vec2mat(LineInfo(i).Cmatrix,LineInfo(i).numPhases)*CapConversionFactor * LineInfo(i).length;
    LineP.(NewLineNames{i})=LineInfo(i).nodes;
    LineZ.(NewLineNames{i})=LineR.(NewLineNames{i}) + 1i*LineX.(NewLineNames{i});
    LineY.(NewLineNames{i})=inv(LineZ.(NewLineNames{i}));
    LineFlow.(NewLineNames{i})=LineInfo(i).bus1PhasePowerReal(LineInfo(i).nodes) + 1i*(LineInfo(i).bus1PhasePowerReactive(LineInfo(i).nodes));
end

Bus=struct();
for i = 1:length(BusInfo)
    Bus(i).vref=slack_voltage;
    Bus(i).phase=(BusInfo(rmap(i)).nodes)';
    Bus(i).kVBase=sqrt(3)*BusInfo(rmap(i)).kVBase;
    Bus(i).vm=nonzeros(BusInfo(rmap(i)).phaseVoltagesPU);
    Bus(i).py=zeros(1,3);
    Bus(i).qy=zeros(1,3);
    Bus(i).sy=zeros(1,3);
    Bus(i).yy=zeros(1,3);
    Bus(i).pd=zeros(1,3);
    Bus(i).qd=zeros(1,3);
    Bus(i).yd=zeros(1,3);
    Bus(i).sd=zeros(1,3);
    Bus(i).Ysh=zeros(3,3);
end

Branch=struct();
for i = 1:length(NewLineNames)
    Branch(i).f=fnew(i);
    Branch(i).t=tnew(i);
    Branch(i).R=LineR.(strcat('l_',string(i)))/zbase;
    Branch(i).X=LineX.(strcat('l_',string(i)))/zbase;
    Branch(i).phase=LineP.(strcat('l_',string(i)));
    Branch(i).Y=LineY.(strcat('l_',string(i)))*zbase;
    Branch(i).Z=LineZ.(strcat('l_',string(i)))/zbase;
    Branch(i).Bc=LineB.(strcat('l_',string(i)))*zbase;
    Branch(i).S=LineFlow.(strcat('l_',string(i)))/(1000*sbase/3);
end

% Adding the branch C to corresponding buses
for br=1:length(Branch)
    temp=zeros(3,3);
    from=Branch(br).f;
    temp(Branch(br).phase',Branch(br).phase')=Branch(br).Bc/2;
    Bus(from).Ysh=Bus(from).Ysh+ temp;
    
    temp=zeros(3,3);
    to=Branch(br).t;
    temp(Branch(br).phase',Branch(br).phase')=Branch(br).Bc/2;
    Bus(to).Ysh=Bus(to).Ysh+ temp;
end
%% Managing the loads 
loads=getLoadInfo(DSSObj);
for i = 1:length(loads)
    loadbusname=loads(i).busName;
    nodes=loads(i).nodes;
    id=nmap(BusNametoId.(strcat('b',loadbusname)));
    pdperphase=loads(i).kW/loads(i).numPhases;
    qdperphase=loads(i).kvar/loads(i).numPhases;
    % Checking the Y connected  Loads 
    if (loads(i).isDelta==0)
        % Constant Power Loads 
        if (strcmpi(loads(i).model,'dssloadconstpq'))
            temp=zeros(1,3);
            temp(nodes)=pdperphase; 
            Bus(id).py = Bus(id).py + temp;
            temp=zeros(1,3);
            temp(nodes)=qdperphase; 
            Bus(id).qy = Bus(id).qy + temp;
        % Constant Impedance Loads 
        elseif (strcmpi(loads(i).model,'dssloadconstz'))
            if (loads(i).numPhases>1)
               zload=(loads(i).kV/sqrt(3))^2*1000/conj(pdperphase+1i*qdperphase)/zbase;
            else
               zload=(loads(i).kV)^2*1000/conj(pdperphase+1i*qdperphase)/zbase; 
            end
            yload=1/zload;
            temp=zeros(1,3);
            temp(nodes)=yload;
            Bus(id).yy=Bus(id).yy + temp;
        else
            error('Load Model Not Supported Right Now.');
        end
    else  % Delta Connected Loads 
       if (strcmpi(loads(i).model,'dssloadconstpq'))
            if (length(nodes) == 3 && loads(i).numPhases ==3)
                temp=zeros(1,3);
                temp(nodes)=pdperphase; 
                Bus(id).pd = Bus(id).pd + temp;
                temp=zeros(1,3);
                temp(nodes)=qdperphase; 
                Bus(id).qd = Bus(id).qd + temp;
            elseif (length(nodes)==2 && loads(i).numPhases ==1)
                temp=zeros(1,3);
                temp(delta_index(nodes(1),nodes(2)))=pdperphase;
                Bus(id).pd = Bus(id).pd + temp;
                temp=zeros(1,3);
                temp(delta_index(nodes(1),nodes(2)))=qdperphase;
                Bus(id).qd = Bus(id).qd + temp;
            end
            Bus(id).sd = Bus(id).pd + 1i* Bus(id).qd;
%             Bus(id).sd = (Baseparam.delta2Yconvmatrix*Bus(id).sd')';
            
       elseif (strcmpi(loads(i).model,'dssloadconstz'))
            if (length(nodes) == 3 && loads(i).numPhases ==3)
                zload=(loads(i).kV)^2*1000/conj(pdperphase+1i*qdperphase)/zbase;
                yload=1/zload;
                temp=zeros(1,3);
                temp(nodes)=yload;
                Bus(id).yd=Bus(id).yd + temp;
            elseif (length(nodes)==2 && loads(i).numPhases ==1)
                zload=(loads(i).kV)^2*1000/conj(pdperphase+1i*qdperphase)/zbase;
                yload=1/zload;
                temp=zeros(1,3);
                temp(delta_index(nodes(1),nodes(2)))=yload;
                Bus(id).yd=Bus(id).yd + temp;               
            else
                warning ('Other loads not supported, loads ignored.');
            end
       else
           error('Load Model Not Supported Right Now.');
       end
    end
   
end

for i = 1:length(BusNames)
    Bus(i).vref=slack_voltage;
    Bus(i).py=Bus(i).py(Bus(i).phase)/(1000*sbase/3);
    Bus(i).qy=Bus(i).qy(Bus(i).phase)/(1000*sbase/3);
    Bus(i).sy=Bus(i).py + 1i* Bus(i).qy ;
%     Bus(i).Ysh = Bus(i).Ysh + diag(Bus(i).yy);
    Bus(i).Ysh=Bus(i).Ysh(Bus(i).phase',Bus(i).phase');
    Bus(i).yy=Bus(i).yy(Bus(i).phase);
    Bus(i).sd=Bus(i).sd/(1000*sbase);
end
% BusVoltage = [Bus.vm];

%% Error Checking 
for i = 1:length(Bus)
    if (length(Bus(i).Ysh) ~= length(Bus(i).phase))
        disp(i)
    end
end

end

