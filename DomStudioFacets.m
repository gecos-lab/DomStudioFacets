function DomStudioFacets(~)
% @ 2020 by Andrea Bistacchi, distributed under the GNU AGPL v3.0 license.
%
% Function used for the analysis of data output from CloudCompare facets 
%
% Last update 2018/09/12

% initialize
clear all
close all
set(0,'DefaultFigureWindowStyle','docked')
rad = pi/180;
deg = 180/pi;

% ______________________________
% 0 - import data and initialize

% load CSV file with facets data
[file, path] = uigetfile('*.csv');
filename = [path file];
[path,file,~] = fileparts(filename);
inDataTable = importdata(filename,',',1);

% copy data to variables
Index = inDataTable.data(:,1);
Center = inDataTable.data(:,5:7);
Length = inDataTable.data(:,11);
Height = inDataTable.data(:,12);
Surface = inDataTable.data(:,10); % this is the real area of the patch with irregular shape
DipAzimuth = inDataTable.data(:,8);
Dip = inDataTable.data(:,9);

clear inDataTable

% Ndata = number of data points
Ndata = length(Index);

% Calculate strike (for rose diagram)
Strike = (DipAzimuth<90).*(DipAzimuth+270) + (DipAzimuth>=90).*(DipAzimuth-90);
symmetricStrike = [Strike; (Strike<=180).*(Strike+180)+(Strike>180).*(Strike-180)]; % used for rose plot

% Calculate poles to planes
Plunge = 90-Dip;
Trend = DipAzimuth+180;
Trend = Trend-(Trend>360)*360;

% Downwards normals oriented as plunge/trend (data in csv table are strange)
Normal = zeros(Ndata,3);
for i = 1:Ndata
    Normal(i,1) = cos(Plunge(i)*rad)*sin(Trend(i)*rad);
    Normal(i,2) = cos(Plunge(i)*rad)*cos(Trend(i)*rad);
    Normal(i,3) = -sin(Plunge(i)*rad);
end

% ________________________________________________________
% 1 - stereoplot, rose diagram, and orientation statistics

% deviation of individual measurements with respect to first measure (arbitrary, just to take one), as cos(angle) from dot product
deltaFirstNormal = Normal*Normal(1,:)';

% reverse sense for normals pointing opposite to first measure
corrNormal = Normal.*(deltaFirstNormal>=0) - Normal.*(deltaFirstNormal<0);

% "bi-directional" Fisher mean, valid for orientation data (just orientation, undefined sense)
sumNormal = sum(corrNormal);
R = sqrt(sumNormal*sumNormal');
meanNormal = sumNormal/R;

% check if meanNormal points downwards, otherwise reverse it and all corrNormal vectors
downward = (meanNormal(3)<0) + (meanNormal(3)==0) - (meanNormal(3)>0);
meanNormal = meanNormal * downward;
corrNormal = corrNormal * downward;

% deviation of each vector (as orientation only) with respect to mean, in degrees
deltaNormal = acos(corrNormal*meanNormal')*deg;

% calculate plunge/trend of meanNormalDouble
meanP = asin(-meanNormal(3))*deg;
meanT = atan2(meanNormal(1),meanNormal(2))*deg;
meanT = meanT+(meanT<0)*360;
meanP = meanP*(meanP>=0)-meanP*(meanP<0);
meanT = meanT-(meanP<0)*180;
meanT = meanT+(meanT<0)*360;
meanT = meanT-(meanT>360)*360;

% calculate dip/azimut of mean plane
meanDip = 90-meanP;
meanDir = meanT+180;
meanDir = meanDir+(meanDir<0)*360;
meanDir = meanDir-(meanDir>360)*360;

%Fisher's K
if Ndata < 16
    fisherK = Ndata/(Ndata-R) * (1 - 1/Ndata)^2; % Fisher, 1953; Tarling, 1971; R.W. Allmendinger, 2016, Stereonet Manual (Davis, 2002 says somebody use (Ndata-2)
else
    fisherK = (Ndata-1)/(Ndata-R); % Fisher, 1953; Tarling, 1971; R.W. Allmendinger, 2016, Stereonet Manual (Davis, 2002 says somebody use (Ndata-2)
end

prob = 0.05; % probability 0.05 indicates that we are 95% confident that the mean vector for a hypothetical large population of measurements would be within a circle of radius alpha95 degrees from the mean vector of our much more limited sample.
alpha95 = deg * acos( 1 - (Ndata-R)/R * ((1/prob)^(1/(Ndata-1)) -1 ) ); % confidence cone half-apical angle (cone axis = mean)

% create contouring grid
yGrid = 1;
gridPlunge(yGrid) = 90;
gridTrend(yGrid) = 0;

for i = 1:10
    m = 6*i;
    radius = i/10;
    DeltaPhi = 360/m;
    for j = 1:m
        phi = j*DeltaPhi;
        yGrid = yGrid+1;
        gridTrend(yGrid) = phi;
        theta = 2 * asin(radius/sqrt(2));
        gridPlunge(yGrid) = 90 - (theta*deg);
    end
end

% direction cosines of node vectors
normalGrid = zeros(331,3);
for i = 1:331
    normalGrid(i,1) = cos(gridPlunge(i)*rad) * sin(gridTrend(i)*rad);
    normalGrid(i,2) = cos(gridPlunge(i)*rad) * cos(gridTrend(i)*rad);
    normalGrid(i,3) = -sin(gridPlunge(i)*rad);
end

% equal area projection coordinates of nodes
for i = 1:331
    xGrid(i) = sqrt(2) * sin((90-gridPlunge(i))*rad/2) .* sin(gridTrend(i)*rad);
    yGrid(i) = sqrt(2) * sin((90-gridPlunge(i))*rad/2) .* cos(gridTrend(i)*rad);
end

%define counting cone
%key = 2.0; cone = (acos(Ndata/(Ndata+key^2)))*deg;
cone = 8.1096144559941786958201832872484; %%% http://en.wikipedia.org/wiki/Solid_angle

% count point density at all 331 nodes
density = zeros(1,331);
parfor i = 1:331
    %dot product of node vectors and data vectors
    for j = 1:Ndata
        deltaTheta = (acos(normalGrid(i,1)*Normal(j,1) + normalGrid(i,2)*Normal(j,2) + normalGrid(i,3)*Normal(j,3)))*deg;
        if deltaTheta <= cone
            density(i) = density(i)+1;
        end
    end
end

% convert density into percent of Ndata
density = density./Ndata*100;

% Kalsbeek contouring

%set parameters:
%number of x and y nodes for grid
%make these smaller if you get an error about excceding the max array size
nx = 50;
ny = nx;

%number of contours:
ci = 5;

%Grid the data linearly:
xStep = abs((min(xGrid)-max(xGrid))/nx);
yStep = abs((min(yGrid)-max(yGrid))/ny);
xGridArray = [min(xGrid):xStep:max(xGrid)];
yGridArray = [min(yGrid):yStep:max(yGrid)];
[xi,yi,zi] = griddata(xGrid,yGrid,density,xGridArray,yGridArray');

v = floor(linspace(min(density),max(density),ci));

% STEREOPLOT
figure(1)

%draw primitive circle
cmap =[ones(1,10)' linspace(1,0,10)' linspace(1,0,10)'];

radius = 1;

% grid
theta = linspace(0,2*pi,180);
x = radius * cos(theta);
y = radius * sin(theta);
plot(x,y,'k-')
hold on

% cosmetics
axis equal
axis off
title({['N. facets = ' num2str(Ndata) ' - Mean = ' num2str(meanDip) '/' num2str(meanDir) ' - Fisher K = ' num2str(fisherK) ] ;...
    ['half-apical angle of the 95% confidence cone = ' num2str(alpha95)];...
    'Concentrations % of total per 1% area' ;...
    ['Maximum concentration = ' num2str(max(density),3) '%']});

% plot density image
imagesc(xGridArray',yGridArray',zi)
colormap(cmap)
plot(x,y,'k-')
plot(0,0,'k+')
plot(radius,0,'k+')
plot(0,radius,'k+')
plot(-radius,0,'k+')
plot(0,-radius,'k+')
plot(0,0,'k+')

% plot density contours
[c,h] = contour(xi,yi,zi,v,'-k');

% label contours
clabel(c,h);

% % print to file
% print('-djpeg',[filename '_contour.jpg'])

% ROSE DIAGRAM
figure(2);

% geologic rose plot inpired by earth_rose.m by Cameron Sparr
% in turn inspired by wind_rose.m
% The code to change the axis labels from a mathematical system to
% N,S,E,W were written by Jordan Rosenthal in a forum post:
% http://www.mathworks.com/matlabcentral/newsreader/author/4601

D = mod(90 - symmetricStrike, 360)*pi/180;

rose(D, 36);   % rose with bins at 10° increment

hHiddenText = findall(gca,'type','text');
Angles = 0 : 30 : 330;
hObjToDelete = zeros( length(Angles)-4, 1 );
k = 0;
for ang = Angles
    hObj = findall(hHiddenText,'string',num2str(ang));
    switch ang
        case 0
            set(hObj,'string','E','HorizontalAlignment','Left');
        case 90
            set(hObj,'string','N','VerticalAlignment','Bottom');
        case 180
            set(hObj,'string','W','HorizontalAlignment','Right');
        case 270
            set(hObj,'string','S','VerticalAlignment','Top');
        otherwise
            k = k + 1;
            hObjToDelete(k) = hObj;
    end
end
delete( hObjToDelete(hObjToDelete~=0) );

title({['N. facets = ' num2str(Ndata) ' - Mean = ' num2str(meanDip) '/' num2str(meanDir) ' - Fisher K = ' num2str(fisherK) ] ;...
    ['half-apical angle of the 95% confidence cone = ' num2str(alpha95)];...
    'Concentrations % of total per 1% area' ;...
    ['Maximum concentration = ' num2str(max(density),3) '%']});


% ___________________
% 2 - plot deviations

% interpolate average outcrop plane with pincipal components
[coeff,~,~] = pca(Center);
outcropNormal = coeff(:,3);

% outcrop size (just used for visualization)
outcropHeight = -1;
while outcropHeight<=0
    outcropHeight = input('outcrop height along dip [100]: ');
    if isempty(outcropHeight), outcropHeight = 100; end
end

outcropLength = -1;
while outcropLength<=0
    outcropLength = input('outcrop length along strike [100]: ');
    if isempty(outcropLength), outcropLength = 100; end
end
disp(' ');

% plunge/trend of outcrop plane
outcropP = asin(-outcropNormal(3))*deg;
outcropT = atan2(outcropNormal(1),outcropNormal(2))*deg;
outcropT = outcropT+(outcropT<0)*360;
outcropP = outcropP*(outcropP>=0)-outcropP*(outcropP<0);
outcropT = outcropT-(outcropP<0)*180;
outcropT = outcropT+(outcropT<0)*360;
outcropT = outcropT-(outcropT>360)*360;
disp(['outcropP = ' num2str(outcropP)]);
disp(['outcropT = ' num2str(outcropT)]);

% dip/dip azimuth of outcrop plane
outcropDip = 90-outcropP;
outcropDipAzimuth = outcropT+180;
outcropDipAzimuth = outcropDipAzimuth+(outcropDipAzimuth<0)*360;
outcropDipAzimuth = outcropDipAzimuth-(outcropDipAzimuth>360)*360;
disp(['outcropDip = ' num2str(outcropDip)]);
disp(['outcropDipAzimuth = ' num2str(outcropDipAzimuth)]);

% outcrop strike
outcropStrike = (outcropDipAzimuth<90).*(outcropDipAzimuth+270) + (outcropDipAzimuth>=90).*(outcropDipAzimuth-90);
disp(['outcropStrike = ' num2str(outcropStrike)]);

% downdip vector with length = 1/2 Height
halfOutcropPlaneDip(1) = cos(outcropDip*rad)*sin(outcropDipAzimuth*rad)*outcropHeight/2;
halfOutcropPlaneDip(2) = cos(outcropDip*rad)*cos(outcropDipAzimuth*rad)*outcropHeight/2;
halfOutcropPlaneDip(3) = sin(outcropDip*rad)*outcropHeight/2;

% alongstrike vector with length = 1/2 Length
halfOutcropPlaneStrike(1) = sin(outcropStrike*rad)*outcropLength/2;
halfOutcropPlaneStrike(2) = cos(outcropStrike*rad)*outcropLength/2;
halfOutcropPlaneStrike(3) = 0;

% four corners of outcrop best fit plane, which is centerd in outcropCenter
outcropCorner1 = - halfOutcropPlaneStrike - halfOutcropPlaneDip;
outcropCorner2 = - halfOutcropPlaneStrike + halfOutcropPlaneDip;
outcropCorner3 = + halfOutcropPlaneStrike + halfOutcropPlaneDip;
outcropCorner4 = + halfOutcropPlaneStrike - halfOutcropPlaneDip;

% intersection (unit vector) of mean fracture plane and outcrop plane
facetOutcropIntersection = cross(meanNormal,outcropNormal);
facetOutcropIntersectionMod = sqrt(facetOutcropIntersection(1)^2+facetOutcropIntersection(2)^2+facetOutcropIntersection(3)^2);
facetOutcropIntersection = facetOutcropIntersection/facetOutcropIntersectionMod;

% plunge/trend of intersection
intersectionPlaneP = asin(-facetOutcropIntersection(3))*deg;
intersectionPlaneT = atan2(facetOutcropIntersection(1),facetOutcropIntersection(2))*deg;
intersectionPlaneT = intersectionPlaneT+(intersectionPlaneT<0)*360;
intersectionPlaneP = intersectionPlaneP*(intersectionPlaneP>=0)-intersectionPlaneP*(intersectionPlaneP<0);
intersectionPlaneT = intersectionPlaneT-(intersectionPlaneP<0)*180;
intersectionPlaneT = intersectionPlaneT+(intersectionPlaneT<0)*360;
intersectionPlaneT = intersectionPlaneT-(intersectionPlaneT>360)*360;
disp(['intersectionPlaneP = ' num2str(intersectionPlaneP)]);
disp(['intersectionPlaneT = ' num2str(intersectionPlaneT)]);

% normal to mean normal and intersection
projectionPlaneNormal = cross(facetOutcropIntersection,meanNormal);
projectionPlaneNormalMod = sqrt(projectionPlaneNormal(1)^2+projectionPlaneNormal(2)^2+projectionPlaneNormal(3)^2);
projectionPlaneNormal = projectionPlaneNormal/projectionPlaneNormalMod;

% plunge/trend of projection plane
projectionPlaneP = asin(-projectionPlaneNormal(3))*deg;
projectionPlaneT = atan2(projectionPlaneNormal(1),projectionPlaneNormal(2))*deg;
projectionPlaneT = projectionPlaneT+(projectionPlaneT<0)*360;
projectionPlaneP = projectionPlaneP*(projectionPlaneP>=0)-projectionPlaneP*(projectionPlaneP<0);
projectionPlaneT = projectionPlaneT-(projectionPlaneP<0)*180;
projectionPlaneT = projectionPlaneT+(projectionPlaneT<0)*360;
projectionPlaneT = projectionPlaneT-(projectionPlaneT>360)*360;
disp(['projectionPlaneP = ' num2str(projectionPlaneP)]);
disp(['projectionPlaneT = ' num2str(projectionPlaneT)]);

% dip/dip azimuth of projection plane
projectionPlaneDip = 90-projectionPlaneP;
projectionPlaneDipAzimuth = projectionPlaneT+180;
projectionPlaneDipAzimuth = projectionPlaneDipAzimuth+(projectionPlaneDipAzimuth<0)*360;
projectionPlaneDipAzimuth = projectionPlaneDipAzimuth-(projectionPlaneDipAzimuth>360)*360;
disp(['projectionPlaneDip = ' num2str(projectionPlaneDip)]);
disp(['projectionPlaneDipAzimuth = ' num2str(projectionPlaneDipAzimuth)]);

% projection plane strike
projectionPlaneStrike = (projectionPlaneDipAzimuth<90).*(projectionPlaneDipAzimuth+270) + (projectionPlaneDipAzimuth>=90).*(projectionPlaneDipAzimuth-90);
disp(['projectionPlaneStrike = ' num2str(projectionPlaneStrike)]);

% downdip vector with length = 1/2 Height
halfProjectionPlaneDip(1) = cos(projectionPlaneDip*rad)*sin(projectionPlaneDipAzimuth*rad)*outcropHeight/2;  % uses outcrop height
halfProjectionPlaneDip(2) = cos(projectionPlaneDip*rad)*cos(projectionPlaneDipAzimuth*rad)*outcropHeight/2;
halfProjectionPlaneDip(3) = sin(projectionPlaneDip*rad)*outcropHeight/2;

% alongstrike vector with length = 1/2 Length
halfProjectionPlaneStrike(1) = sin(projectionPlaneStrike*rad)*outcropLength/2;  % uses outcrop length
halfProjectionPlaneStrike(2) = cos(projectionPlaneStrike*rad)*outcropLength/2;
halfProjectionPlaneStrike(3) = 0;

% four corners of outcrop best fit plane, which is centerd in outcropCenter
projectionPlaneCorner1 = - halfProjectionPlaneStrike - halfProjectionPlaneDip;
projectionPlaneCorner2 = - halfProjectionPlaneStrike + halfProjectionPlaneDip;
projectionPlaneCorner3 = + halfProjectionPlaneStrike + halfProjectionPlaneDip;
projectionPlaneCorner4 = + halfProjectionPlaneStrike - halfProjectionPlaneDip;

% for each FACET, downdip vector with length = 1/2 Height
halfDip = zeros(Ndata,3);
for i = 1:Ndata
    halfDip(i,1) = cos(Dip(i)*rad)*sin(DipAzimuth(i)*rad)*Height(i)/2;
    halfDip(i,2) = cos(Dip(i)*rad)*cos(DipAzimuth(i)*rad)*Height(i)/2;
    halfDip(i,3) = -sin(Dip(i)*rad)*Height(i)/2;
end

% for each FACET, alongstrike vector with length = 1/2 Length
halfStrike = zeros(Ndata,3);
for i = 1:Ndata
    halfStrike(i,1) = sin(Strike(i)*rad)*Length(i)/2;
    halfStrike(i,2) = cos(Strike(i)*rad)*Length(i)/2;
    halfStrike(i,3) = 0;
end

% for each FACET, four corners of each fracture in reference frame with center = outcropCenter
corner1 = Center - halfStrike - halfDip - outcropCenter;
corner2 = Center - halfStrike + halfDip - outcropCenter;
corner3 = Center + halfStrike + halfDip - outcropCenter;
corner4 = Center + halfStrike - halfDip - outcropCenter;

% plot fracture rectangles in 3D
figure(4); hold on; axis equal

plot3([corner1(:,1) corner2(:,1) corner3(:,1) corner4(:,1) corner1(:,1)]',...
    [corner1(:,2) corner2(:,2) corner3(:,2) corner4(:,2) corner1(:,2)]',...
    [corner1(:,3) corner2(:,3) corner3(:,3) corner4(:,3) corner1(:,3)]');
quiver3(0,0,0,meanNormal(1)*outcropLength/3,meanNormal(2)*outcropLength/3,meanNormal(3)*outcropLength/3,'LineWidth',2,'Color','cyan');

plot3([outcropCorner1(1) outcropCorner2(1) outcropCorner3(1) outcropCorner4(1) outcropCorner1(1)]',...
    [outcropCorner1(2) outcropCorner2(2) outcropCorner3(2) outcropCorner4(2) outcropCorner1(2)]',...
    [outcropCorner1(3) outcropCorner2(3) outcropCorner3(3) outcropCorner4(3) outcropCorner1(3)]',...
    'LineWidth',2,'Color','magenta');
quiver3(0,0,0,outcropNormal(1)*outcropLength/3,outcropNormal(2)*outcropLength/3,outcropNormal(3)*outcropLength/3,'LineWidth',2,'Color','magenta');

plot3([projectionPlaneCorner1(1) projectionPlaneCorner2(1) projectionPlaneCorner3(1) projectionPlaneCorner4(1) projectionPlaneCorner1(1)]',...
    [projectionPlaneCorner1(2) projectionPlaneCorner2(2) projectionPlaneCorner3(2) projectionPlaneCorner4(2) projectionPlaneCorner1(2)]',...
    [projectionPlaneCorner1(3) projectionPlaneCorner2(3) projectionPlaneCorner3(3) projectionPlaneCorner4(3) projectionPlaneCorner1(3)]',...
    'LineWidth',2,'Color','green');
quiver3(0,0,0,projectionPlaneNormal(1)*outcropLength/3,projectionPlaneNormal(2)*outcropLength/3,projectionPlaneNormal(3)*outcropLength/3,'LineWidth',2,'Color','green');

quiver3(0,0,0,facetOutcropIntersection(1)*outcropLength/3,facetOutcropIntersection(2)*outcropLength/3,facetOutcropIntersection(3)*outcropLength/3,'LineWidth',2,'Color','red');

title({'3D view of fracture facets - mean normal in cyan, outcrop in magenta,';...
       'projection plane in green, intersection of outcrop and mean facet plane in red'})
xlabel('East'); ylabel('North'); zlabel('Z');

% ______________________________________________
% 5 - project facets onto projection plane in 3D

% find projections of corner points of fractures
corner1Projected = zeros(Ndata,3);
corner2Projected = zeros(Ndata,3);
corner3Projected = zeros(Ndata,3);
corner4Projected = zeros(Ndata,3);

parfor i = 1:Ndata
    corner1Projected(i,:) = corner1(i,:) + (-corner1(i,:) * projectionPlaneNormal') * projectionPlaneNormal;
    corner2Projected(i,:) = corner2(i,:) + (-corner2(i,:) * projectionPlaneNormal') * projectionPlaneNormal;
    corner3Projected(i,:) = corner3(i,:) + (-corner3(i,:) * projectionPlaneNormal') * projectionPlaneNormal;
    corner4Projected(i,:) = corner4(i,:) + (-corner4(i,:) * projectionPlaneNormal') * projectionPlaneNormal;
end

% plot projected fracture rectangles in 3D
figure(5); hold on; axis equal

plot3([corner1Projected(:,1) corner2Projected(:,1) corner3Projected(:,1) corner4Projected(:,1) corner1Projected(:,1)]',...
    [corner1Projected(:,2) corner2Projected(:,2) corner3Projected(:,2) corner4Projected(:,2) corner1Projected(:,2)]',...
    [corner1Projected(:,3) corner2Projected(:,3) corner3Projected(:,3) corner4Projected(:,3) corner1Projected(:,3)]');
quiver3(0,0,0,meanNormal(1)*outcropLength/3,meanNormal(2)*outcropLength/3,meanNormal(3)*outcropLength/3,'LineWidth',2,'Color','cyan');

plot3([outcropCorner1(1) outcropCorner2(1) outcropCorner3(1) outcropCorner4(1) outcropCorner1(1)]',...
    [outcropCorner1(2) outcropCorner2(2) outcropCorner3(2) outcropCorner4(2) outcropCorner1(2)]',...
    [outcropCorner1(3) outcropCorner2(3) outcropCorner3(3) outcropCorner4(3) outcropCorner1(3)]',...
    'LineWidth',2,'Color','magenta');
quiver3(0,0,0,outcropNormal(1)*outcropLength/3,outcropNormal(2)*outcropLength/3,outcropNormal(3)*outcropLength/3,'LineWidth',2,'Color','magenta');

plot3([projectionPlaneCorner1(1) projectionPlaneCorner2(1) projectionPlaneCorner3(1) projectionPlaneCorner4(1) projectionPlaneCorner1(1)]',...
    [projectionPlaneCorner1(2) projectionPlaneCorner2(2) projectionPlaneCorner3(2) projectionPlaneCorner4(2) projectionPlaneCorner1(2)]',...
    [projectionPlaneCorner1(3) projectionPlaneCorner2(3) projectionPlaneCorner3(3) projectionPlaneCorner4(3) projectionPlaneCorner1(3)]',...
    'LineWidth',2,'Color','green');
quiver3(0,0,0,projectionPlaneNormal(1)*outcropLength/3,projectionPlaneNormal(2)*outcropLength/3,projectionPlaneNormal(3)*outcropLength/3,'LineWidth',2,'Color','green');
quiver3(0,0,0,facetOutcropIntersection(1)*outcropLength/3,facetOutcropIntersection(2)*outcropLength/3,facetOutcropIntersection(3)*outcropLength/3,'LineWidth',2,'Color','red');

title({'3D view of projected fracture facets - mean normal in cyan, outcrop in magenta,';...
       'projection plane in green, intersection of outcrop and mean facet plane in red'})
xlabel('East'); ylabel('North'); zlabel('Z');

% __________________________
% 6 - facets traces 3D -> 2D

UVcorner1 = zeros(Ndata,2);
UVcorner2 = zeros(Ndata,2);
UVcorner3 = zeros(Ndata,2);
UVcorner4 = zeros(Ndata,2);

% U coorddinate along scanline unit vector
parfor i = 1:Ndata
    UVcorner1(i,1) = dot(corner1Projected(i,:),meanNormal);
    UVcorner2(i,1) = dot(corner2Projected(i,:),meanNormal);
    UVcorner3(i,1) = dot(corner3Projected(i,:),meanNormal);
    UVcorner4(i,1) = dot(corner4Projected(i,:),meanNormal);
end

% V coordinate along intersection unit vector
parfor i = 1:Ndata
    UVcorner1(i,2) = dot(corner1Projected(i,:),facetOutcropIntersection);
    UVcorner2(i,2) = dot(corner2Projected(i,:),facetOutcropIntersection);
    UVcorner3(i,2) = dot(corner3Projected(i,:),facetOutcropIntersection);
    UVcorner4(i,2) = dot(corner4Projected(i,:),facetOutcropIntersection);
end

% plot projected fracture rectangles in 3D
figure(6); hold on; axis equal; set(gca,'XDir','reverse'); set(gca,'YDir','reverse')

plot([UVcorner1(:,1) UVcorner2(:,1) UVcorner3(:,1) UVcorner4(:,1) UVcorner1(:,1)]',...
    [UVcorner1(:,2) UVcorner2(:,2) UVcorner3(:,2) UVcorner4(:,2) UVcorner1(:,2)]');

title('2D view of projected fracture facets')
xlabel('mean normal axis (U)'); ylabel('outcrop - mean facet intersection (V)');

% ______________________________________________________________
% 7 - simplify facet traces to segments with just two end-points

UVbottom = zeros(Ndata,2);

for i = 1:Ndata
    Uarray = [UVcorner1(i,1) UVcorner2(i,1) UVcorner3(i,1) UVcorner4(i,1)];
    Varray = [UVcorner1(i,2) UVcorner2(i,2) UVcorner3(i,2) UVcorner4(i,2)];
    
    UVtop(i,1) = Uarray(Varray==max(Varray));
    UVtop(i,2) = Varray(Varray==max(Varray));
    UVbottom(i,1) = Uarray(Varray==min(Varray));
    UVbottom(i,2) = Varray(Varray==min(Varray));
end

% plot projected fracture rectangles in 3D
figure(7); hold on; axis equal; set(gca,'XDir','reverse'); set(gca,'YDir','reverse')

plot([UVtop(:,1) UVbottom(:,1)]',...
    [UVtop(:,2) UVbottom(:,2)]',"LineWidth",2);

title('2D view of projected fracture traces')
xlabel('mean normal axis (U)'); ylabel('outcrop - mean facet intersection (V)');

% _________________________________________________________________________
% 8 - input number of scanlines and find intersections for each facet trace

% input number of scanlines
disp(' ');
Nscan = -1;
while Nscan<=0
    Nscan = input('input number of scanlines [100]: ');
    if isempty(Nscan), Nscan = 100; end
    Nscan = round(Nscan);
end
disp(' ');

% scanline max and min U are given by dataset max and min U plus some tolerance
maxU = max([UVtop(:,1)' UVbottom(:,1)']);
minU = min([UVtop(:,1)' UVbottom(:,1)']);
scanMaxU = maxU + outcropLength*0.005;
scanMinU = minU - outcropLength*0.005;

% scanline step along V is given dividing maxV - minV minus some tolerance by number of scanlines - 1
maxV = max(UVtop(:,2));
minV = min(UVbottom(:,2));
scanMaxV = maxV - outcropHeight*0.005;
scanMinV = minV + outcropHeight*0.005;
scanStep = (scanMaxV - scanMinV)/(Nscan-1);

intersectionCount = zeros(1,Nscan);
spacing = [];

for i = 1:Nscan
    scanV = scanMinV + scanStep*(i-1);
    plot([scanMaxU scanMinU],[scanV scanV],"Color",[0.5 0.5 0.5]);
    recordedIntersection = [];
    
    for j = 1:Ndata
        [intersectionU,intersectionV] = polyxpoly([UVtop(j,1) UVbottom(j,1)],[UVtop(j,2) UVbottom(j,2)],[scanMaxU scanMinU],[scanV scanV]);
        if isfinite(intersectionU)
            plot(intersectionU,intersectionV,'k.');
            intersectionCount(i) = intersectionCount(i)+1;
            recordedIntersection = [recordedIntersection intersectionU];
        end
    end
    
    if intersectionCount(i)>1
        disp(['Scanline ' num2str(i) ' - ' num2str(intersectionCount(i)) ' intersections']);
        recordedIntersection = sort(recordedIntersection);
        thisScanSpacing = recordedIntersection(2:end) - recordedIntersection(1:end-1);
        %disp(num2str(thisScanSpacing));
        spacing = [spacing thisScanSpacing];
    end
end

spacingPrctile = prctile(spacing,0:5:100);
spacingMean = mean(spacing);
spacingMedian = median(spacing);
spacingMode = mode(spacing);
spacingStDev = std(spacing);

figure(8);
histogram(spacing,'BinWidth',1,'Normalization','pdf');
grid on 
xlabel('Spacing'); ylabel('Frequency');
title({['Spacing histogram from ' num2str(Nscan) ' scanlines'],...
       ['Mean = ' num2str(spacingMean) '  Mode = ' num2str(spacingMode) '  St Dev = ' num2str(spacingStDev)],...
       ['Percentiles:'],...
       ['  0% = ' num2str(spacingPrctile(1),3) '   5% = ' num2str(spacingPrctile(2),3) '  10% = ' num2str(spacingPrctile(3),3) '  15% = ' num2str(spacingPrctile(4),3) '  20% = ' num2str(spacingPrctile(5),3)],...
       [' 25% = ' num2str(spacingPrctile(6),3) '  30% = ' num2str(spacingPrctile(7),3) '  35% = ' num2str(spacingPrctile(8),3) '  40% = ' num2str(spacingPrctile(9),3) '  45% = ' num2str(spacingPrctile(10),3)],...
       [' 50% = ' num2str(spacingPrctile(11),3) '  55% = ' num2str(spacingPrctile(12),3) '  60% = ' num2str(spacingPrctile(13),3) '  65% = ' num2str(spacingPrctile(14),3) '  70% = ' num2str(spacingPrctile(15),3)],...
       [' 75% = ' num2str(spacingPrctile(16),3) '  80% = ' num2str(spacingPrctile(17),3) '  85% = ' num2str(spacingPrctile(18),3) '  90% = ' num2str(spacingPrctile(19),3) '  95% = ' num2str(spacingPrctile(20),3)],...
       ['100% = ' num2str(spacingPrctile(21),3) ]});

figure(9);
histogram(spacing,'BinWidth',1,'Normalization','cdf');
xlabel('Spacing'); ylabel('Cumulative Frequency');
title({['Cumulative frequency from ' num2str(Nscan) ' scanlines'],...
       ['Mean = ' num2str(spacingMean) '  Mode = ' num2str(spacingMode) '  St Dev = ' num2str(spacingStDev)],...
       ['Percentiles:'],...
       ['  0% = ' num2str(spacingPrctile(1),3) '   5% = ' num2str(spacingPrctile(2),3) '  10% = ' num2str(spacingPrctile(3),3) '  15% = ' num2str(spacingPrctile(4),3) '  20% = ' num2str(spacingPrctile(5),3)],...
       [' 25% = ' num2str(spacingPrctile(6),3) '  30% = ' num2str(spacingPrctile(7),3) '  35% = ' num2str(spacingPrctile(8),3) '  40% = ' num2str(spacingPrctile(9),3) '  45% = ' num2str(spacingPrctile(10),3)],...
       [' 50% = ' num2str(spacingPrctile(11),3) '  55% = ' num2str(spacingPrctile(12),3) '  60% = ' num2str(spacingPrctile(13),3) '  65% = ' num2str(spacingPrctile(14),3) '  70% = ' num2str(spacingPrctile(15),3)],...
       [' 75% = ' num2str(spacingPrctile(16),3) '  80% = ' num2str(spacingPrctile(17),3) '  85% = ' num2str(spacingPrctile(18),3) '  90% = ' num2str(spacingPrctile(19),3) '  95% = ' num2str(spacingPrctile(20),3)],...
       ['100% = ' num2str(spacingPrctile(21),3) ]});

% ______________________________
% 9 - save figures for reporting
for i=1:9
    figure(i);
    %savefig([path '\' file '_fig_' num2str(i) '.fig']);
    saveas(gcf,[path '\' file '_fig_' num2str(i) '.jpg'])
end

% save spacing
spacing = spacing';
save([path '\' file '_spacing.txt'],'spacing','-ascii');

disp('Files saved');

