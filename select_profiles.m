function [float_ids, float_profs] = select_profiles(lon_lim,lat_lim,...
    start_date,end_date,varargin)
% select_profiles  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   [float_ids, float_profs] = select_profiles(lon_lim,lat_lim,...
%       start_date,end_date,varargin)
%
% DESCRIPTION:
%   This function returns the indices of profiles and floats that match
%   the given criteria (spatial, temporal, sensor availability).
%   It does not download any files.
%
% INPUTS:
% lon_lim    : longitude limits
% lat_lim    : latitude limits
%            * Latitude and longitude limits can be input as either
%            two element vectors ([LON1 LON2], [LAT1 LAT2]) for maximum
%            and minimum limits, or as same-sized vectors with at least
%            3 elements for vertices of a polygon
%            * Longitude can be input in either the -180 to 180 degrees
%            format or 0 to 360 degrees format (or even in any other
%            360 degree range that encloses all the desired longitude
%            values, e.g., [-20 200])
%            * Either or both values can be '[]' to indicate the full range
% start_date : start date
% end_date   : end date
%            * Dates should be in one of the following formats:
%            [YYYY MM DD HH MM SS] or [YYYY MM DD]
%            * Both values can be '[]' to indicate the full range
%
% OPITONAL INPUTS (key,value pairs):
% 'outside', 'none' 'time' 'space' 'both': By default, only float profiles
%           that are within both the temporal and spatial constraints are
%           returned ('none'); specify to also maintain profiles outside
%           the temporal constraints ('time'), spatial constraints
%           ('space'), or both constraints ('both')
% 'sensor', 'SENSOR_TYPE': By default, all floats within the lon/lat/time
%           limits are considered. This option allows the selection by 
%           sensor type. Available are: PRES, PSAL, TEMP, DOXY, BBP,
%           BBP470, BBP532, BBP700, TURBIDITY, CP, CP660, CHLA, CDOM,
%           NITRATE, BISULFIDE, PH_IN_SITU_TOTAL, DOWN_IRRADIANCE,
%           DOWN_IRRADIANCE380, DOWN_IRRADIANCE412, DOWN_IRRADIANCE443, 
%           DOWN_IRRADIANCE490, DOWN_IRRADIANCE555, DOWN_IRRADIANCE670, 
%           UP_RADIANCE, UP_RADIANCE412, UP_RADIANCE443, UP_RADIANCE490,
%           UP_RADIANCE555, DOWNWELLING_PAR, DOXY2, DOXY3
%           (Currently, only one sensor type can be selected.)
% 'ocean', ocean: Valid choices are 'A' (Atlantic), 'P' (Pacific), and
%           'I' (Indian). This selection is in addition to the specified
%           longitude and latitude limits. (To select all floats and 
%           profiles from one ocean basin, leave lon_lim and lat_lim
%           empty.)
% 'mode',mode: Valid modes are 'R' (real-time), 'A' (adjusted), and
%           'D', in any combination. Only profiles with the selected
%           mode(s) will be listed in float_profs.
%           Default is 'RAD' (all modes).
%           If 'sensor' option is not used, the 'mode' option is ignored.
%
% OUTPUTS:
%   float_ids   : array with the WMO IDs of all matching floats
%   float_profs : cell array with the per-float indices of all matching profiles
%
% AUTHORS: 
%   J. Sharp, H. Frenzel, A. Fassbender (NOAA-PMEL), N. Buzby (UW),
%   J. Plant, T. Maurer, Y. Takeshita (MBARI), D. Nicholson (WHOI),
%   and A. Gray (UW)
%
% CITATION:
%   H. Frenzel*, J. Sharp*, A. Fassbender, N. Buzby, J. Plant, T. Maurer,
%   Y. Takeshita, D. Nicholson, A. Gray, 2021. BGC-Argo-Mat: A MATLAB
%   toolbox for accessing and visualizing Biogeochemical Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.4971318.
%   (*These authors contributed equally to the code.)
%
% LICENSE: bgc_argo_mat_license.m
%
% DATE: DECEMBER 1, 2021  (Version 1.1)

global Float Settings Sprof;

% make sure Settings is initialized
if isempty(Settings)
    initialize_argo();
end

% set defaults
outside = 'none'; % if set, removes profiles outside time/space constraints
sensor = []; % default: use all available profiles
ocean = []; % default: any ocean basin
mode = 'RAD';

% parse optional arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'outside')
        outside = varargin{i+1};
    elseif strcmpi(varargin{i}, 'sensor')
        sensor = varargin{i+1};
    elseif strcmpi(varargin{i}, 'ocean')
        ocean = upper(varargin{i+1}(1));
    elseif strcmpi(varargin{i}, 'mode')
        mode = varargin{i+1};
    end
end

% check if specified sensor exists
if ~isempty(sensor)
    if ~any(strcmp(sensor, Settings.avail_vars))
        warning('unknown sensor: %s', sensor)
        pause(3)
        sensor = [];
    end
end

% check if specified ocean is correct
if ~isempty(ocean) && ~contains('API', ocean)
    warning('no such ocean: %s', ocean)
    ocean = [];
end

% check if specified data modes are correct
new_mode= '';
for i = 1:length(mode)
    if contains('RAD', mode(i))
        new_mode = strcat(new_mode, mode(i));
    else
        warning('no such mode: %s', mode(i))
    end
end
if isempty(new_mode)
    mode = 'ADR';
else
    mode = sort(new_mode); % standard order enables strcmp later
end

% make sure Sprof is initialized
if isempty(Sprof)
    initialize_argo();
end

% fill in the blanks if needed
if isempty(lon_lim)
    lon_lim = [-180, 180];
end
if isempty(lat_lim)
    lat_lim = [-90, 90];
end
if isempty(start_date)
    start_date = [1995, 1, 1];
end
if isempty(end_date)
    end_date = [2038, 1, 19];
end

% ADJUST INPUT DATES TO DATENUM FORMAT
dn1 = datenum(start_date); 
dn2 = datenum(end_date);

% GET INDEX OF PROFILES WITHIN USER-SPECIFIED GEOGRAPHIC POLYGON
inpoly = get_inpolygon(Sprof.lon,Sprof.lat,lon_lim,lat_lim);
if isempty(inpoly) || ~any(inpoly)
    warning('no matching profiles found')
    float_ids = [];
    float_profs = [];
    return
end

% Find index of dates that are within the time window
date_inpoly = datenum(Sprof.date(inpoly), 'yyyymmddHHMMSS');
indate_poly = date_inpoly >= dn1 & date_inpoly <= dn2;
% now create an indate array of 0s/1s that has the same 
% size as inpoly so that it can be used in the & operations below
indate = zeros(size(inpoly));
all_floats = 1:length(inpoly);
sel_floats_space = all_floats(inpoly);
indate(sel_floats_space(indate_poly)) = 1;

% SELECT BY SENSOR
if isempty(sensor)
    has_sensor = ones(size(indate)); % no sensor was selected
else
    has_sensor = contains(Sprof.sens, sensor);
end
if ~any(has_sensor)
    warning('no data found for sensor %s', sensor);
end

% select by ocean basin
if isempty(ocean)
    is_ocean = ones(size(indate)); % no ocean was selected
else
    is_ocean = strcmp(Sprof.ocean, ocean);
end

% select by data mode
if isempty(mode) || strcmp(mode, 'ADR')
    has_mode = ones(size(inpoly));
else
    has_mode = zeros(size(inpoly));
    is_good = inpoly & indate & has_sensor & is_ocean;
    idx = all_floats(is_good);
    for i = 1:length(idx)
        pos = strcmp(split(Sprof.sens(idx(i))), sensor);
        if any(pos)
            str = cell2mat(Sprof.data_mode(idx(i)));
            has_mode(idx(i)) = contains(mode, str(pos));
        end
    end
end

% perform selection
all_prof = 1:length(indate);
profiles = all_prof(inpoly & indate & has_sensor & is_ocean & has_mode);
float_ids = unique(Sprof.wmo(profiles));

% download Sprof files if necessary
good_float_ids = download_multi_floats(float_ids);
 
% the information from the index file is only used for an initial
% filtering of floats, the actual information from the Sprof files
% is used in a second step
float_ids = good_float_ids;
float_profs = cell(length(good_float_ids), 1);

for fl = 1:length(good_float_ids)
    filename = sprintf('%s%d_Sprof.nc', Settings.prof_dir, ...
        good_float_ids(fl));
    [n_prof, n_param] = get_dims(filename);
    fl_idx = find(Float.wmoid == good_float_ids(fl), 1);
    n_prof_exp = Float.prof_idx2(fl_idx) - Float.prof_idx1(fl_idx) + 1;
    if n_prof_exp > n_prof
        warning(['The index file lists %d profiles for float %d, ', ...
            'but the Sprof file has only %d profiles.'], ...
            n_prof_exp, good_float_ids(fl), n_prof)
    end
    lon = ncread(filename, 'LONGITUDE');
    lat = ncread(filename, 'LATITUDE');
    juld = ncread(filename, 'JULD');
    if ~isempty(sensor) && ~strcmp(mode, 'ADR')
        params = ncread(filename, 'PARAMETER');
        param_names = cell(n_param, 1);
        for p = 1:n_param
            param_names{p} = strtrim(params(:,p,1,1)');
        end
        param_idx = find(strcmp(param_names, sensor), 1);
        data_mode = ncread(filename, 'PARAMETER_DATA_MODE');
    end
    date = datenum(juld) + datenum([1950 1 1]);

    inpoly = get_inpolygon(lon,lat,lon_lim,lat_lim);
    indate = date >= dn1 & date <= dn2;

    if isempty(sensor)
        has_sensor = ones(size(inpoly));
    else
        param = ncread(filename, 'PARAMETER');
        has_sensor = zeros(size(inpoly));
        for p = 1:n_prof
           has_sensor(p) = any(strcmp(cellstr(param(:,:,1,p)'), sensor));
        end
    end
    if isempty(ocean)
        is_ocean = ones(size(inpoly));
    else
        sprof_loc = Sprof.lon + 1i * Sprof.lat;
        this_loc = lon + 1i * lat;
        [~,idx] = min(abs(bsxfun(@minus,sprof_loc(:),this_loc(:).')));
        is_ocean = strcmp(Sprof.ocean(idx),ocean);
        is_ocean(isnan(this_loc)) = 0;
    end
    if isempty(mode) || strcmp(mode, 'ADR')
        has_mode = ones(size(inpoly));
    else
        has_mode = zeros(size(inpoly));
        for m = 1:length(mode)
            has_mode = has_mode + (data_mode(param_idx,:)' == mode(m));
        end
    end
    % now apply the given constraints
    all_prof = 1:length(inpoly);
    if strcmp(outside, 'none')
        float_profs{fl} = all_prof(inpoly & indate & has_sensor & ...
            is_ocean & has_mode);
    elseif strcmp(outside, 'time') % must meet space constraint
        float_profs{fl} = all_prof(inpoly & has_sensor & is_ocean & has_mode);
    elseif strcmp(outside, 'space') % must meet time constraint
        float_profs{fl} = all_prof(indate & has_sensor & is_ocean & has_mode);
    elseif strcmp(outside, 'both') % no time or space constraint
        float_profs{fl} = all_prof(has_sensor & is_ocean & has_mode);
    else
        warning('no such setting for "outside": %s', outside)
        float_profs{fl} = [];
    end
 
    if isempty(float_profs{fl})
        warning('no matching profiles found for float %d', good_float_ids(fl))
        float_ids(float_ids == good_float_ids(fl)) = [];
    end
end

float_profs(cellfun(@isempty, float_profs)) = [];
