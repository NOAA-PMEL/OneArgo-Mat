function float_ids = select_profiles_per_type(Profiles, ...
    lon_lim, lat_lim, dn1, dn2, interp_ll, sensor, ocean)
% select_profiles_per_type  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   float_ids = select_profiles_per_type(Profiles, ...
%       lon_lim, lat_lim, dn1, dn2, interp_ll, sensor, ocean)
%
% DESCRIPTION:
%   This function returns the IDs of floats that match
%   the given criteria (spatial, temporal, sensor availability). It is a
%   helper function for select_profiles and does not select by all
%   criteria specified by the user with arguments to select_profiles.
%   Note: due to inconsistencies between index and Profiles files, selecting
%   by data mode  is not performed in the first round, only in the second
%   round (in select_profiles)
%
% INPUTS:
%   lon_lim : longitude limits
%   lat_lim : latitude limits
%            * Latitude and longitude limits can be input as either
%            two element vectors ([LON1 LON2], [LAT1 LAT2]) for maximum
%            and minimum limits, or as same-sized vectors with at least
%            3 elements for vertices of a polygon
%            * Longitude can be input in either the -180 to 180 degrees
%            format or 0 to 360 degrees format (or even in any other
%            360 degree range that encloses all the desired longitude
%            values, e.g., [-20 200])
%   dn1 : start date (in datenum format)
%   dn2 : end date (in datenum format)
%   interp_ll : if 'yes', missing lon/lat values will be interpolated
%   sensor : This option allows the selection by
%           sensor type. Available are: PRES, PSAL, TEMP, DOXY, BBP,
%           BBP470, BBP532, BBP700, TURBIDITY, CP, CP660, CHLA, CDOM,
%           NITRATE, BISULFIDE, PH_IN_SITU_TOTAL, DOWN_IRRADIANCE,
%           DOWN_IRRADIANCE380, DOWN_IRRADIANCE412, DOWN_IRRADIANCE443,
%           DOWN_IRRADIANCE490, DOWN_IRRADIANCE555, DOWN_IRRADIANCE670,
%           UP_RADIANCE, UP_RADIANCE412, UP_RADIANCE443, UP_RADIANCE490,
%           UP_RADIANCE555, DOWNWELLING_PAR, CNDC, DOXY2, DOXY3, BBP700_2
%           (Full list can be displayed with the list_sensors function.)
%           Multiple sensors can be entered as a cell array, e.g.:
%           {'DOXY';'NITRATE'}
%   ocean : Valid choices are 'A' (Atlantic), 'P' (Pacific), and
%           'I' (Indian). This selection is in addition to the specified
%           longitude and latitude limits.
%
% OUTPUTS:
%   float_ids : array with the WMO IDs of all matching floats
%
% AUTHORS:
%   H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL), N. Buzby (UW)
%
% CITATION:
%   H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2022. OneArgo-Mat:
%   A MATLAB toolbox for accessing and visualizing Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.6588042
%
% LICENSE: oneargo_mat_license.m
%
% DATE: JUNE 1, 2022  (Version 1.0.1)

float_ids = []; % default return value

% GET INDEX OF PROFILES WITHIN USER-SPECIFIED GEOGRAPHIC POLYGON
inpoly = get_inpolygon(Profiles.lon,Profiles.lat,lon_lim,lat_lim);

% if interpolation of missing values is requested, include all missing
% positions for now; actual interpolation will be performed in the
% second round of matching after reading the Profiles files
if strncmpi(interp_ll, 'yes', 1)
    inpoly(isnan(Profiles.lon)) = 1;
end

if isempty(inpoly) || ~any(inpoly)
    warning('no matching profiles found')
    return
end

% Find index of dates that are within the time window
date_inpoly = datenum(Profiles.date(inpoly), 'yyyymmddHHMMSS');
indate_poly = date_inpoly >= dn1 & date_inpoly <= dn2;
% now create an indate array of 0s/1s that has the same
% size as inpoly so that it can be used in the & operations below
indate = zeros(size(inpoly));
all_floats = 1:length(inpoly);
sel_floats_space = all_floats(inpoly);
indate(sel_floats_space(indate_poly)) = 1;

% select by sensor
has_sensor = ones(size(indate));
if ~isempty(sensor)
    for i = 1:length(sensor)
        has_sensor = has_sensor & cellfun(@(x) ...
            any(strcmp(x, sensor{i})), Profiles.split_sens);
    end
end
if ~any(has_sensor)
    warning('no profiles found that have all specified sensors')
    return
end

% select by ocean basin
if isempty(ocean)
    is_ocean = ones(size(indate)); % no ocean was selected
else
    is_ocean = strcmp(Profiles.ocean, ocean);
end

% perform selection
all_prof = 1:length(indate);
profiles = all_prof(inpoly & indate & has_sensor & is_ocean);
float_ids = unique(Profiles.wmo(profiles));
