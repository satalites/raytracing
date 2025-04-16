% Read building data from OSM file
filename = "ucla.osm";
buildings = readgeotable(filename, Layer="buildingparts");
basemapName = "osm";
url = "a.tile.openstreetmap.org/${z}/${x}/${y}.png"; 
addCustomBasemap(basemapName, url);

% Tx location - lat and lon corresponding to latitude and longitude)
txLat = 34.072604;
txLon = -118.442273;

% Building ROI - region where buildings will be loaded in
buildingLatROI = [34.07491; 34.06384; 34.06384; 34.07491];
buildingLonROI = [-118.44850; -118.44850; -118.43940; -118.43940];
buildingLatROI(end+1) = buildingLatROI(1);
buildingLonROI(end+1) = buildingLonROI(1);

% Defining region where users will be placed
userLatROI = [34.072506; 34.071933; 34.071933; 34.072506];
userLonROI = [-118.444891; -118.444891; -118.440787; -118.440787];

% Spacing between users (m)
spacingMeters = 0.5;

ROI = geopolyshape(buildingLatROI, buildingLonROI);

[latmin, latmax] = bounds(buildingLatROI);
[lonmin, lonmax] = bounds(buildingLonROI);

shape = buildings.Shape;
clipped = geoclip(shape, [latmin, latmax], [lonmin, lonmax]);
idxInsideROI = clipped.NumRegions > 0;
buildingsROI = buildings(idxInsideROI,:);

% set building properties
for row = 1:height(buildingsROI)
    buildingsROI.Material(row) = "brick";
    buildingsROI.Color(row) = "#AA4A44";
end

% hide fig(can be changed by setting Visible=true)
viewer = siteviewer(Buildings=buildingsROI, Basemap=basemapName);

tx = txsite(Latitude=txLat, Longitude=txLon, TransmitterFrequency=28e9);

pm = propagationModel("raytracing", MaxNumReflections=5);

[userLatMin, userLatMax] = bounds(userLatROI);
[userLonMin, userLonMax] = bounds(userLonROI);

% Convert spacing to deg (approx)
latSpacing = spacingMeters / 111000;  % 1 degree latitude ~ 111 km
lonSpacing = spacingMeters / (111000 * cosd(mean(userLatROI)));

% Determine # of grid pts
latRange = userLatMax - userLatMin;
lonRange = userLonMax - userLonMin;
latPoints = max(1, floor(latRange / latSpacing) + 1);
lonPoints = max(1, floor(lonRange / lonSpacing) + 1);

latGrid = linspace(userLatMin, userLatMin + (latPoints-1)*latSpacing, latPoints);
lonGrid = linspace(userLonMin, userLonMin + (lonPoints-1)*lonSpacing, lonPoints);
[lonMesh, latMesh] = meshgrid(lonGrid, latGrid);
userPositions = [latMesh(:), lonMesh(:)];
numUsers = size(userPositions, 1);

fprintf('\nStarting ray tracing analysis for %d users with spacing of %.1f meters...\n', numUsers, spacingMeters);
rayProperties = struct();

% for userIdx = 1:numUsers
%     userLat = userPositions(userIdx, 1);
%     userLon = userPositions(userIdx, 2);
%     rxUser = rxsite(Latitude=userLat, Longitude=userLon);
    
%     % perform ray trace
%     rays = raytrace(tx, rxUser, pm);
%     userRays = rays{1};
    
%     % sort by path loss
%     if ~isempty(userRays)
%         [~, sortedIdx] = sort([userRays.PathLoss], 'ascend');
%         topRayIndices = sortedIdx(1:min(5, numel(userRays)));
%     else
%         topRayIndices = [];
%     end
    
%     % Calculate grid coordinates (row, col)
%     [~, rowIdx] = min(abs(latGrid - userLat));
%     [~, colIdx] = min(abs(lonGrid - userLon));
    
%     % store properties for each ray
%     userField = sprintf('User_%d', userIdx);
%     rayProperties.(userField) = struct('RayIdx', {}, 'AoD_Az', {}, 'AoD_El', {}, ...
%         'AoA_Az', {}, 'AoA_El', {}, 'PhaseShift', {}, ...
%         'Position', {}, 'PathLoss', {}, 'GridCoord', {});
    
%     for k = 1:length(topRayIndices)
%         rayIdx = topRayIndices(k);
%         ray = userRays(rayIdx);
%         AoD = ray.AngleOfDeparture;
%         AoA = ray.AngleOfArrival;
        
%         rayProps = struct('RayIdx', rayIdx, ...
%              'AoD_Az', AoD(1), ...
%              'AoD_El', AoD(2), ...
%              'AoA_Az', AoA(1), ...
%              'AoA_El', AoA(2), ...
%              'PhaseShift', ray.PhaseShift, ...
%              'Position', [userLat, userLon], ...
%              'PathLoss', ray.PathLoss, ...
%              'GridCoord', [rowIdx, colIdx]);
%         rayProperties.(userField)(k) = rayProps;
%     end
    
%     fprintf('Completed ray tracing for User %d of %d\n', userIdx, numUsers);
% end

% % save data
% save('spacedGridRayProperties.mat', '-struct', 'rayProperties');
% save('userGridInfo.mat', 'userPositions', 'latPoints', 'lonPoints', 'spacingMeters', ...
%     'latGrid', 'lonGrid', 'tx');

% clear viewer
% set(0, 'DefaultFigureVisible', 'on');
% fprintf('\nAll operations completed.\n');