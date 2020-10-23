%% plotGoogleMap
% Plots a google map on the current axes using the Google Static Maps API
%
%% Syntax
%  h = plotGoogleMap(Property, Value,...);
%  [lonVec latVec imag] = plotGoogleMap(Property, Value,...);
%
%% Description
% Plots the google map on the current axes given the input properties
% selected
%
%% Inputs
% * *|Property|* - property name from the list below along with the
% * *|Value|* for the property.  The default for each porperty is in
% parenthesis.
% * -- *'MapType'*         - ('roadmap')  Type of map to return. Any of [roadmap, 
%                                  satellite, terrain, hybrid) See the Google Maps API for
%                                  more information. 
% * -- *'Alpha'* (1)         - (0-1) Transparency level of the map (0 is fully
%                                  transparent). While the map is always
%                                  moved to the bottom of the plot (i.e. will
%                                  not hide previously drawn items), this can
%                                  be useful in order to increase readability
%                                  if many colors are plotted (using SCATTER
%                                  for example).
% * -- *'Marker'*            - The marker argument is a text string with fields
%                                  conforming to the Google Maps API. The
%                                  following are valid examples:
%                                  '43.0738740,-70.713993' (dflt midsize orange marker)
%                                  '43.0738740,-70.713993,blue' (midsize blue marker)
%                                  '43.0738740,-70.713993,yellowa' (midsize yellow
%                                  marker with label "A")
%                                  '43.0738740,-70.713993,tinyredb' (tiny red marker
%                                  with label "B")
%
%% Outputs
% * *|h|* - Handle to the plotted map
% * *|lonVect|* - Vector of Longidute coordinates (WGS84) of the image
% * *|latVect|* - Vector of Latidute coordinates (WGS84) of the image
% * *|imag|* - Image matrix (height,width,3) of the map
%
%% References:
% http://www.mathworks.com/matlabcentral/fileexchange/24113
% http://www.maptiler.org/google-maps-coordinates-tile-bounds-projection/
% http://developers.google.com/maps/documentation/staticmaps/
%
%  Acknowledgement to Val Schmidt for his submission of get_google_map.m
%  Acknowledgement to Zohar Bar-Yehuda for his submission of plot_google_map.mp
%
%% Copyright
% Copyright (c) 2010, Zohar Bar-Yehuda
% Copyright (c) 2010, Val Schmidt
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
% 
%% Example
% Plot a map showing some capitals in Europe:
%%
% lat = [48.8708   51.5188   41.9260   40.4312   52.523   37.982];
% lon = [2.4131    -0.1300    12.4951   -3.6788    13.415   23.715];
% plot(lon,lat,'.r','MarkerSize',20)
% plotGoogleMap
%

function varargout = plotGoogleMap(varargin)


% Handle input arguments
height = 640;
width = 640;
scale = 2;
maptype = 'hybrid';
alphaData = 1;
hold on

if nargin >= 2
    for idx = 1:2:length(varargin)
        switch lower(varargin{idx})
            case 'maptype'
                maptype = varargin{idx+1};
            case 'alpha'
                alphaData = varargin{idx+1};
            otherwise
                error(['Unrecognized variable: ' varargin{idx}])
        end
    end
end

curAxis = axis;
if isequal(curAxis,[0 1 0 1]) % probably an empty figure
    % display US map
    curAxis = [-130 -65 25 50];
    axis(curAxis)
end

% Store paramters in axis handle (for auto refresh callbacks)
ud = get(gca, 'UserData');
ud.gmap_params = varargin;
set(gca, 'UserData', ud);

% adjust current axis limit to avoid strectched maps
[xExtent,yExtent] = latLonToMeters(curAxis(3:4), curAxis(1:2) );
xExtent = diff(xExtent); % just the size of the span
yExtent = diff(yExtent);
if xExtent > yExtent
    % enlarge the X extent
    centerY = mean(curAxis(3:4));
    spanY = (curAxis(4)-curAxis(3))/2;
    curAxis(3) = centerY-spanY*xExtent/yExtent;
    curAxis(4) = centerY+spanY*xExtent/yExtent;
elseif yExtent > xExtent
    % enlarge the Y extent
    centerX = mean(curAxis(1:2));
    spanX = (curAxis(2)-curAxis(1))/2;
    curAxis(1) = centerX-spanX*yExtent/xExtent;
    curAxis(2) = centerX+spanX*yExtent/xExtent;
end
axis(curAxis) % update axis as quickly as possible, before downloading new image
drawnow

% Delete previous map from plot (if exists)
if nargout <= 1 % only if in plotting mode
    curChildren = get(gca,'children');
    delete(findobj(curChildren,'tag','gmap'))
end

% Enforce Latitude constraints of EPSG:900913 
if curAxis(3) < -85
    curAxis(3) = -85;
end
if curAxis(4) > 85
    curAxis(4) = 85;
end

% Calculate zoom level for current axis limits
[xExtent,yExtent] = latLonToMeters(curAxis(3:4), curAxis(1:2) );
minResX = diff(xExtent) / width;
minResY = diff(yExtent) / height;
minRes = max([minResX minResY])/2;
tileSize = 256;
initialResolution = 2 * pi * 6378137 / tileSize; % 156543.03392804062 for tileSize 256 pixels
zoomlevel = floor(log2(initialResolution/minRes));

% Enforce valid zoom levels
if zoomlevel < 0 
    zoomlevel = 0;
end
if zoomlevel > 19 
    zoomlevel = 19;
end

% Calculate center coordinate in WGS1984 of four images
lat = (curAxis(4)-curAxis(3))/2+curAxis(3);
lon = (curAxis(2)-curAxis(1))/2+curAxis(1);
lat1 = (curAxis(4)-curAxis(3))/4+curAxis(3);
lon1 = (curAxis(2)-curAxis(1))/4+curAxis(1);
lat2 = (curAxis(4)-curAxis(3))*3/4+curAxis(3);
lon2 = (curAxis(2)-curAxis(1))/4+curAxis(1);
lat3 = (curAxis(4)-curAxis(3))/4+curAxis(3);
lon3 = (curAxis(2)-curAxis(1))*3/4+curAxis(1);
lat4 = (curAxis(4)-curAxis(3))*3/4+curAxis(3);
lon4 = (curAxis(2)-curAxis(1))*3/4+curAxis(1);

% CONSTRUCT QUERY URL
preamble = 'http://maps.googleapis.com/maps/api/staticmap';
location1 = ['?center=' num2str(lat1,10) ',' num2str(lon1,10)];
location2 = ['?center=' num2str(lat2,10) ',' num2str(lon2,10)];
location3 = ['?center=' num2str(lat3,10) ',' num2str(lon3,10)];
location4 = ['?center=' num2str(lat4,10) ',' num2str(lon4,10)];
zoomStr = ['&zoom=' num2str(zoomlevel)];
sizeStr = ['&scale=' num2str(scale) '&size=' num2str(width) 'x' num2str(height)];
maptypeStr = ['&maptype=' maptype ];

filename1 = [cd,'\tmp1.jpg'];
filename2 = [cd,'\tmp2.jpg'];
filename3 = [cd,'\tmp3.jpg'];
filename4 = [cd,'\tmp4.jpg'];
format = '&format=jpg';

url1 = [preamble location1 zoomStr sizeStr maptypeStr format '&sensor=false'];
url2 = [preamble location2 zoomStr sizeStr maptypeStr format '&sensor=false'];
url3 = [preamble location3 zoomStr sizeStr maptypeStr format '&sensor=false'];
url4 = [preamble location4 zoomStr sizeStr maptypeStr format '&sensor=false'];


% Get the image
try
    urlwrite(url1,filename1);
    urlwrite(url2,filename2);
    urlwrite(url3,filename3);
    urlwrite(url4,filename4);
catch % error downloading map
    warning('Unable to download map form Google Servers.\nPossible reasons: no network connection, or quota exceeded (1000 map requests per day).')
    return
end
%pause(1)
try
    [M1 Mcolor1] = imread(filename1);
    [M2 Mcolor2] = imread(filename2);
    [M3 Mcolor3] = imread(filename3);
    [M4 Mcolor4] = imread(filename4);
catch
    pause(1);
    [M1 Mcolor1] = imread(filename1);
    [M2 Mcolor2] = imread(filename2);
    [M3 Mcolor3] = imread(filename3);
    [M4 Mcolor4] = imread(filename4);
end
M1 = cast(M1,'double');
M2 = cast(M2,'double');
M3 = cast(M3,'double');
M4 = cast(M4,'double');
% pause(1)
try
    delete(filename1); % delete temp file
    delete(filename2); % delete temp file
    delete(filename3); % delete temp file
    delete(filename4); % delete temp file
catch
    %oh well
end

width = size(M1,2);
height = size(M1,1);

% Calculate a meshgrid of pixel coordinates in EPSG:900913
centerPixelY = round(height/2);
centerPixelX = round(width/2);
[centerX1,centerY1] = latLonToMeters(lat1, lon1 ); % center coordinates in EPSG:900913
[centerX2,centerY2] = latLonToMeters(lat2, lon2 ); % center coordinates in EPSG:900913
[centerX3,centerY3] = latLonToMeters(lat3, lon3 ); % center coordinates in EPSG:900913
[centerX4,centerY4] = latLonToMeters(lat4, lon4 ); % center coordinates in EPSG:900913

curResolution = initialResolution / 2^zoomlevel/scale; % meters/pixel (EPSG:900913)
xVec1 = centerX1 + ((1:width)-centerPixelX) * curResolution; % x vector
yVec1 = centerY1 + ((height:-1:1)-centerPixelY) * curResolution; % y vector
[xMesh1,yMesh1] = meshgrid(xVec1,yVec1); % construct meshgrid 
xVec2 = centerX2 + ((1:width)-centerPixelX) * curResolution; % x vector
yVec2 = centerY2 + ((height:-1:1)-centerPixelY) * curResolution; % y vector
[xMesh2,yMesh2] = meshgrid(xVec2,yVec2); % construct meshgrid 
xVec3 = centerX3 + ((1:width)-centerPixelX) * curResolution; % x vector
yVec3 = centerY3 + ((height:-1:1)-centerPixelY) * curResolution; % y vector
[xMesh3,yMesh3] = meshgrid(xVec3,yVec3); % construct meshgrid 
xVec4 = centerX4 + ((1:width)-centerPixelX) * curResolution; % x vector
yVec4 = centerY4 + ((height:-1:1)-centerPixelY) * curResolution; % y vector
[xMesh4,yMesh4] = meshgrid(xVec4,yVec4); % construct meshgrid 

% convert meshgrid to WGS1984
[lonMesh1,latMesh1] = metersToLatLon(xMesh1,yMesh1);
[lonMesh2,latMesh2] = metersToLatLon(xMesh2,yMesh2);
[lonMesh3,latMesh3] = metersToLatLon(xMesh3,yMesh3);
[lonMesh4,latMesh4] = metersToLatLon(xMesh4,yMesh4);

imag1 = M1/255;
imag2 = M2/255;
imag3 = M3/255;
imag4 = M4/255;

imag1 = imag1(latMesh1(:,1)<lat,lonMesh1(1,:)<lon,:);
lonMesh1f = lonMesh1(latMesh1(:,1)<lat,lonMesh1(1,:)<lon);
latMesh1f = latMesh1(latMesh1(:,1)<lat,lonMesh1(1,:)<lon);
imag2 = imag2(latMesh2(:,1)>lat,lonMesh2(1,:)<lon,:);
lonMesh2f = lonMesh2(latMesh2(:,1)>lat,lonMesh2(1,:)<lon);
latMesh2f = latMesh2(latMesh2(:,1)>lat,lonMesh2(1,:)<lon);
imag3 = imag3(latMesh3(:,1)<lat,lonMesh3(1,:)>lon,:);
lonMesh3f = lonMesh3(latMesh3(:,1)<lat,lonMesh3(1,:)>lon);
latMesh3f = latMesh3(latMesh3(:,1)<lat,lonMesh3(1,:)>lon);
imag4 = imag4(latMesh4(:,1)>lat,lonMesh4(1,:)>lon,:);
lonMesh4f = lonMesh4(latMesh4(:,1)>lat,lonMesh4(1,:)>lon);
latMesh4f = latMesh4(latMesh4(:,1)>lat,lonMesh4(1,:)>lon);
minDim = min([size(imag1,1),size(imag2,1),size(imag3,1),size(imag4,1),size(imag1,2),size(imag2,2),size(imag3,2),size(imag4,2)]);
imag = [imag2,imag4;imag1,imag3];
lonMesh = [lonMesh2f, lonMesh4f; lonMesh1f, lonMesh3f];
latMesh = [latMesh2f, latMesh4f; latMesh1f, latMesh3f];

if nargout <= 1 % plot map

    % display image
    h = image(lonMesh(1,:),latMesh(:,1),1-(1-imag)*alphaData);
    set(gca,'YDir','normal')
    
    %h = surface(lonMesh,latMesh,zeros(size(latMesh)),1-(1-imag)*alphaData,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct','DisplayName','map');
    
    set(h,'tag','gmap')
    uistack(h,'bottom') % move map to bottom (so it doesn't hide previously drawn annotations)
    axis(curAxis) % restore original zoom
    if nargout == 1
        varargout{1} = h;
    end
    
    % if auto-refresh mode - override zoom callback to allow autumatic 
    % refresh of map upon zoom actions.
    figHandle = gca;
    while ~strcmpi(get(figHandle, 'Type'), 'figure')
        % Recursively search for parent figure in case axes are in a panel
        figHandle = get(figHandle, 'Parent');
    end
    
    zoomHandle = zoom(gca);   
    panHandle = pan(figHandle); % This isn't ideal, doesn't work for contained axis    
    set(zoomHandle,'ActionPostCallback',@update_google_map);
    set(panHandle, 'ActionPostCallback', @update_google_map);

else % don't plot, only return map
    varargout{1} = lonMesh;
    varargout{2} = latMesh;
    varargout{3} = imag;
end

function update_google_map(obj,evd)
% callback function for auto-refresh
drawnow;
try
    axHandle = evd.Axes;
catch ex
    % Event doesn't contain the correct axes. Panic!
    axHandle = gca;
end
ud = get(axHandle, 'UserData');
if isfield(ud, 'gmap_params')
    params = ud.gmap_params;
    plotGoogleMap(params{:});
end

% Coordinate transformation functions

function [lon,lat] = metersToLatLon(x,y)
% Converts XY point from Spherical Mercator EPSG:900913 to lat/lon in WGS84 Datum
originShift = 2 * pi * 6378137 / 2.0; % 20037508.342789244
lon = (x ./ originShift) * 180;
lat = (y ./ originShift) * 180;
lat = 180 / pi * (2 * atan( exp( lat * pi / 180)) - pi / 2);

function [x,y] = latLonToMeters(lat, lon )
% Converts given lat/lon in WGS84 Datum to XY in Spherical Mercator EPSG:900913"
originShift = 2 * pi * 6378137 / 2.0; % 20037508.342789244
x = lon * originShift / 180;
y = log(tan((90 + lat) * pi / 360 )) / (pi / 180);
y = y * originShift / 180;
    