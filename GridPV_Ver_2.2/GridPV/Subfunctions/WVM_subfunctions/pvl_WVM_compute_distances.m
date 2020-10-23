function [dist]=pvl_WVM_compute_distances(lonin,latin,type,MW,PVdensity)

%types = 
%'square' (square PV plant with specified PV density and MWs';
    %latin and lonin are single values
%'polygon'(polygon plant shape with veriticies specified); 
    %latin and lonin are vectors of the polygon verticies 
%'discrete' (e.g., to replicate aggregate output of multiple point sensors)
    %latin and lonin are vectors of the discrete locations
%MW = plant MWs   
%PVdensity = W PV installed per m2 across the plant (only needed for

%% determine amount of memory available
userview=memory;
maxbytes=userview.MaxPossibleArrayBytes;
maxentries=maxbytes./8;

%% Convert Lat/Lon to x-y with distances in meters

R=6371008.7714; %mean radius of Earth, in meters, from "Geodetic Reference System" (1980) by Moritz (http://download.springer.com/static/pdf/409/art%253A10.1007%252Fs001900050278.pdf?auth66=1416866651_cf2d69f45ebfc42880b8986b0f14de40&ext=.pdf).
meters_per_degree_lat=R*pi/180;
meters_per_degree_lon=R*cosd(mean(latin))*pi/180;

X=meters_per_degree_lon.*lonin;
Y=meters_per_degree_lat.*latin;


%% Square Plant
if strcmp(type,'square')==1
    %compute length of plant
    length1=round(sqrt(MW*10^6./PVdensity)); %compute length in meters across of square plant in x-direction
    
    %check to see if computer has enough memory to handle all points at 1m
    %spacing
    maxlength_xsquare=sqrt(maxentries);
    if length1<maxlength_xsquare/10 %if less than 10% of max allowed length, set spacing = 1m 
        spacing=1; %use 1 meter spacing 
        [xsquare,ysquare]=meshgrid(X-length1/2:spacing:X+length1/2,Y-length1./2:spacing:Y+length1./2);
    else
        spacing=round(length1./(maxlength_xsquare/10));
        warning(['Spacing between simulation points was adjusted to ' num2str(spacing) 'm (instead of 1m) due to lack of memory. This may lead to smoother output than intended.'])
        [xsquare,ysquare]=meshgrid(X-length1/2:spacing:X+length1/2,Y-length1./2:spacing:Y+length1./2);
    end
    
    xinds=reshape(xsquare,1,[]);
    yinds=reshape(ysquare,1,[]);
    
end  
    
%% Polygon Plant
if strcmp(type,'polygon')
parea=polyarea(X,Y);
    PVdensity=MW*10^6./parea;
    fprintf(['The PV density based on the polygon inputs was ' sprintf('%2.0f',PVdensity) 'W/m2 \n' ]);
    v2=[min(min(X))-100 max(max(X))+100 min(min(Y))-100 max(max(Y))+100]; %draw square around polygon with some buffer
    
    length2=(max(max(X))+100-(min(min(X))-100)).*(max(max(Y))+100-(min(min(Y))-100));
    
    if length2<maxentries/10
    spacing=1;
    [xgrid,ygrid]=meshgrid(v2(1):spacing:v2(2),v2(3):spacing:v2(4));
    else
       spacing=round(length2./(maxentries./10));
       warning(['Spacing between simulation points was adjusted to ' num2str(spacing) 'm (instead of 1m) due to lack of memory. This may lead to smoother output than intended.'])
       [xgrid,ygrid]=meshgrid(v2(1):spacing:v2(2),v2(3):spacing:v2(4));
    end
    IN=inpolygon(xgrid,ygrid,X,Y);
    xinds=reshape(xgrid,1,[]);
    yinds=reshape(ygrid,1,[]);
    IN2=reshape(IN,1,[]);
    xinds(IN2==0)=[];
    yinds(IN2==0)=[];
end
    
%% Discrete Plant
if strcmp(type,'discrete')==1
    xinds=X;
    yinds=Y;
end
    
%% take a random sampling of xinds since pdist is limited by avaiable memory
%use 1/10 of available memory

if strcmp(type,'discrete')==1
    try
        px=pdist(xinds); %use pdist (Statistics Toolbox) if available
        py=pdist(yinds);
    catch
        px=pdist_vector(xinds); %else, use pdist_vector (slightly slower)
        py=pdist_vector(yinds);
    end
else
    %determine max length
    fn=@(x)abs(x.^2-x-2*maxentries);
    [maxlength_xinds]=fminsearch(fn,3000);
    length_random=min(3000,floor(maxlength_xinds./10)); %take 3000 random samples, or use machine-specified max (10% of available memory)
    
    if length(xinds) > length_random
        if length_random<500
            warning(['The available memory is limited and may lead to inaccurate results. For computing distances between modules, a random sample of all module locations is taken (to reduce memory requirements). Less than 500 samples were taken due to memory constraints. Samples taken: ' num2str(length_random) '.']);
        end
        
        r=randi(length(xinds),length_random,1);
        try
            px=pdist(xinds(r)'); %use pdist (Statistics Toolbox) if available
            py=pdist(yinds(r)'); 
        catch
            px=pdist_vector(xinds(r)'); %else, use pdist_vector (slightly slower)
            py=pdist_vector(yinds(r)');
        end
    else
        try
            px=pdist(xinds'); %use pdist (Statistics Toolbox) if available
            py=pdist(yinds');
        catch
            px=pdist_vector(xinds'); %else, use pdist_vector (slightly slower)
            py=pdist_vector(yinds');
        end
    end
end

dist=sqrt(px.^2+py.^2);    