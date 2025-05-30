% Read building data from OSM file
filename = "ucla.osm";
buildings = readgeotable(filename, Layer="buildingparts");
basemapName = "osm";
url = "a.tile.openstreetmap.org/${z}/${x}/${y}.png"; 
addCustomBasemap(basemapName, url);

% Tx location - lat and lon corresponding to latitude and longitude
txLat = 34.072604;
txLon = -118.442273;

% Building ROI - region where buildings will be loaded in
buildingLatROI = [34.07491; 34.06384; 34.06384; 34.07491];
buildingLonROI = [-118.44850; -118.44850; -118.43940; -118.43940];
% so now we're going to append the first entry in the array to the end of
% it. cait note: i am less familiar with matlab but this does seem to be a
% way of closing a loop of coordinates so that line segments will connect
% properly. cool!
buildingLatROI(end+1) = buildingLatROI(1);
buildingLonROI(end+1) = buildingLonROI(1);

% Defining region where users will be placed
userLatROI = [34.072506; 34.071933; 34.071933; 34.072506];
userLonROI = [-118.444891; -118.444891; -118.440787; -118.440787];
% cait note: VERY similar to the buiding ROI. why the difference?

% Spacing between users (m)
spacingMeters = 10;

% this is what will create our polygon! 
ROI = geopolyshape(buildingLatROI, buildingLonROI);

% cait note: bounding box? we use later for var clipped...
[latmin, latmax] = bounds(buildingLatROI);
[lonmin, lonmax] = bounds(buildingLonROI);

% cait note: wow! let's get into the details, line by line!
% starting with var shape, the geospatial table "buildings," which has been
% read from the .osm file, we extract a subset of the table ("Shape")
shape = buildings.Shape;
% and now, clipped (defined by the bounding box)
% will apply all of our shapes extracted from the buildings variable.
% geoclip itself "clips" the extracted shapes to the box!
clipped = geoclip(shape, [latmin, latmax], [lonmin, lonmax]);
% this creates a mask of which shapes remain after we clip (as in, the
% indices of shapes existing inside clipped area get marked by this mask)
idxInsideROI = clipped.NumRegions > 0;
% now, we can look only at buildings existing in our ROI!
buildingsROI = buildings(idxInsideROI,:);

% set building properties
for row = 1:height(buildingsROI)
    buildingsROI.Material(row) = "brick";
    buildingsROI.Color(row) = "#AA4A44";
end

% hide fig(can be changed by setting Visible=true)
viewer = siteviewer(Buildings=buildingsROI, Basemap=basemapName);

% now we must set a location for our base station
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

% beginning at the user's minimum lat or lon, generate however many points
% specified by latPoints 
% cait note: can we not set the end point as just userLatMax or userLonMax?
latGrid = linspace(userLatMin, userLatMin + (latPoints-1)*latSpacing, latPoints);
lonGrid = linspace(userLonMin, userLonMin + (lonPoints-1)*lonSpacing, lonPoints);

% now, we want to take those points and create a grid
[lonMesh, latMesh] = meshgrid(lonGrid, latGrid);
userPositions = [latMesh(:), lonMesh(:)];
numUsers = size(userPositions, 1);

fprintf('\nStarting ray tracing analysis for %d users with spacing of %.1f meters...\n', numUsers, spacingMeters);
rayProperties = struct();

% now we will iterate through all users and begin ray tracing!

for userIdx = 1:numUsers
    userLat = userPositions(userIdx, 1);
    userLon = userPositions(userIdx, 2);
    % earlier, we set a location for our transmitter. now, we collect the
    % lat/lon position of our user and set our receive location!
    rxUser = rxsite(Latitude=userLat, Longitude=userLon);

    % perform ray trace
    rays = raytrace(tx, rxUser, pm);
    % extract the propagation paths from the cell array
    userRays = rays{1};

    % let's take a peek at what's happening !!
    % show(tx)
    % show(rxUser)
    % if ~isempty(userRays)
    %     plot(userRays)
    % end

    % from here on out, we're going to set up data saving!
    % sort by path loss
    if ~isempty(userRays)
        [~, sortedIdx] = sort([userRays.PathLoss], 'ascend');
        topRayIndices = sortedIdx(1:min(5, numel(userRays)));
    else
        topRayIndices = [];
    end

    % Calculate grid coordinates (row, col)
    [~, rowIdx] = min(abs(latGrid - userLat));
    [~, colIdx] = min(abs(lonGrid - userLon));

    % store properties for each ray
    userField = sprintf('User_%d', userIdx);
    rayProperties.(userField) = struct('RayIdx', {}, 'AoD_Az', {}, 'AoD_El', {}, ...
        'AoA_Az', {}, 'AoA_El', {}, 'PhaseShift', {}, ...
        'Position', {}, 'PathLoss', {}, 'GridCoord', {});

    for k = 1:length(topRayIndices)
        rayIdx = topRayIndices(k);
        ray = userRays(rayIdx);
        AoD = ray.AngleOfDeparture;
        AoA = ray.AngleOfArrival;

        rayProps = struct('RayIdx', rayIdx, ...
             'AoD_Az', AoD(1), ...
             'AoD_El', AoD(2), ...
             'AoA_Az', AoA(1), ...
             'AoA_El', AoA(2), ...
             'PhaseShift', ray.PhaseShift, ...
             'Position', [userLat, userLon], ...
             'PathLoss', ray.PathLoss, ...
             'GridCoord', [rowIdx, colIdx]);
        rayProperties.(userField)(k) = rayProps;
    end

    fprintf('Completed ray tracing for User %d of %d\n', userIdx, numUsers);
end

% save data
save('spacedGridRayProperties.mat', '-struct', 'rayProperties');
save('userGridInfo.mat', 'userPositions', 'latPoints', 'lonPoints', 'spacingMeters', ...
    'latGrid', 'lonGrid', 'tx');

clear viewer
set(0, 'DefaultFigureVisible', 'on');
fprintf('\nAll operations completed.\n');