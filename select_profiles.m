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
%            * Either or both values can be '[]' to indicate the full range
%   start_date : start date
%   end_date : end date
%            * Dates should be in one of the following formats:
%            [YYYY MM DD HH MM SS] or [YYYY MM DD]
%            * Both values can be '[]' to indicate the full range
%
% OPTIONAL INPUTS (key,value pairs):
%   'dac',dac: Select by Data Assimilation Center reponsible for the floats.
%           A single DAC can be entered as a string (e.g.: 'aoml'),
%           multiple DACs can be entered as a cell array (e.g.:
%           {'meds';'incois'}.
%           Valid values are any of: {'aoml'; 'bodc'; 'coriolis'; ...
%           'csio'; 'csiro'; 'incois'; 'jma'; 'kma'; 'kordi'; 'meds'}
%   'depth',depth: Select profiles that reach at least this depth
%           (positive downwards; in db)
%   'floats',floats: Select profiles only from these floats that must
%           match all other criteria
%   'min_num_prof',num_prof: Select only floats that have at least
%           num_prof profiles that meet all other criteria
%   'mode',mode: Valid modes are 'R' (real-time), 'A' (adjusted), and
%           'D', in any combination. Only profiles with the selected
%           mode(s) will be listed in float_profs.
%           Default is 'RAD' (all modes).
%           If multiple sensors are specified, all of them must be in 
%           the selected mode(s).
%           If 'sensor' option is not used, the 'mode' option is ignored.
%   'ocean', ocean: Valid choices are 'A' (Atlantic), 'P' (Pacific), and
%           'I' (Indian). This selection is in addition to the specified
%           longitude and latitude limits. (To select all floats and 
%           profiles from one ocean basin, leave lon_lim and lat_lim
%           empty.)
%   'outside', 'none' 'time' 'space' 'both': By default, only float profiles
%           that are within both the temporal and spatial constraints are
%           returned ('none'); specify to also maintain profiles outside
%           the temporal constraints ('time'), spatial constraints
%           ('space'), or both constraints ('both')
%   'sensor', SENSOR_TYPE: This option allows the selection by 
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
float_ids = [];
float_profs = [];
outside = 'none'; % if set, removes profiles outside time/space constraints
sensor = []; % default: use all profiles that match other criteria
ocean = []; % default: any ocean basin
mode = 'RAD';
dac = [];
floats = [];
depth = [];
min_num_prof = 0;

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
    elseif strcmpi(varargin{i}, 'dac')
        dac = varargin{i+1};
    elseif strcmpi(varargin{i}, 'floats')
        floats = varargin{i+1};
    elseif strcmpi(varargin{i}, 'depth')
        depth = varargin{i+1};
    elseif strcmpi(varargin{i}, 'min_num_prof')
        min_num_prof = varargin{i+1};
    else
        warning('unknown option: %s', varargin{i});
    end
end

% convert requested sensor(s) to cell array if necessary and
% discard unknown sensors
sensor = check_variables(sensor, 'warning', ...
    'unknown sensor will be ignored');

% only use mode if sensor was specified
if isempty(sensor)
    mode = [];
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

% check if specified dac(s) are correct
if ischar(dac)
    dac = cellstr(dac);
end
bad = zeros(length(dac), 1);
for i = 1:length(dac)
    if ~any(strcmp(dac{i}, Settings.dacs))
        warning('no such dac: %s', dac{i});
        bad(i) = 1;
    end
end
dac(bad == 1) = [];

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
has_sensor = ones(size(indate));
if ~isempty(sensor)
    sprof_sens = cellfun(@split, Sprof.sens, 'UniformOutput', false);
    for i = 1:length(sensor)
        has_sensor = has_sensor & cellfun(@(x) ...
            any(strcmp(x, sensor{i})), sprof_sens);
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
    is_ocean = strcmp(Sprof.ocean, ocean);
end

% select by data mode
% Note: due to inconsistencies between index and Sprof files, this
% matching is not performed in the first round, only in the second
% round below (based on Sprof files)
%if strcmp(mode, 'ADR')
has_mode = ones(size(inpoly));
% else
%     has_mode = zeros(size(inpoly));
%     is_good = inpoly & indate & has_sensor & is_ocean;
%     idx = all_floats(is_good);
%     for i = 1:length(idx)
%         pos = strcmp(split(Sprof.sens(idx(i))), sensor);
%         if any(pos)
%             str = cell2mat(Sprof.data_mode(idx(i)));
%             has_mode(idx(i)) = contains(mode, str(pos));
%         end
%     end
% end

% perform selection
all_prof = 1:length(indate);
profiles = all_prof(inpoly & indate & has_sensor & is_ocean & has_mode);
float_ids = unique(Sprof.wmo(profiles));

% check for selected DACs if applicable (DACs are stored by float,
% not by profile)
if ~isempty(dac)
    idx = arrayfun(@(x) find(Float.wmoid==x, 1), float_ids);
    found_dacs = Float.dac(idx);
    uses_dac = ismember(found_dacs, dac);
    float_ids = float_ids(uses_dac);
end

if ~isempty(floats)
    % select only those floats found so far that also were specified
    float_ids = intersect(float_ids, floats);
end

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
    if isempty(depth)
        has_press = ones(size(lon));
    else
        press = ncread(filename, 'PRES');
        has_press = (max(press) >= depth)';
    end
    if ~isempty(sensor) && ~strcmp(mode, 'ADR')
        params = ncread(filename, 'PARAMETER');
        param_names = cell(n_param, 1);
        % find the index of a profile that has the most sensors available
        tmp = sum(sum(params));
        pidx = find(max(tmp) == tmp, 1);
        for p = 1:n_param
            param_names{p} = strtrim(params(:,p,1,pidx)');
        end
        param_idx = zeros(length(sensor), 1);
        for s = 1:length(sensor)
            param_idx(s) = find(strcmp(param_names, sensor{s}), 1);
        end
        data_mode = ncread(filename, 'PARAMETER_DATA_MODE');
    end
    date = datenum(juld) + datenum([1950 1 1]);

    inpoly = get_inpolygon(lon,lat,lon_lim,lat_lim);
    indate = date >= dn1 & date <= dn2;

    has_sensor = ones(size(inpoly));
    if ~isempty(sensor)
        param = ncread(filename, 'PARAMETER');
        for p = 1:n_prof
           for s = 1:length(sensor)
               has_sensor(p) = has_sensor(p) & ...
                   any(strcmp(cellstr(param(:,:,1,p)'), sensor{s}));
           end
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
    if strcmp(mode, 'ADR')
        has_mode = ones(size(inpoly));
    else
        sens_has_mode = zeros(size(inpoly));
        for s = 1:length(sensor)
            for m = 1:length(mode)
                sens_has_mode = sens_has_mode + ...
                    (data_mode(param_idx(s),:)' == mode(m));
            end
        end
        % all sensors must be in one of the specified modes
        has_mode = (sens_has_mode == length(sensor));
    end
    % now apply the given constraints
    all_prof = 1:length(inpoly);
    if strcmp(outside, 'none')
        float_profs{fl} = all_prof(inpoly & indate & has_sensor & ...
            is_ocean & has_mode & has_press);
    elseif strcmp(outside, 'time') % must meet space constraint
        float_profs{fl} = all_prof(inpoly & has_sensor & is_ocean & ...
            has_mode & has_press);
    elseif strcmp(outside, 'space') % must meet time constraint
        float_profs{fl} = all_prof(indate & has_sensor & is_ocean & ...
            has_mode & has_press);
    elseif strcmp(outside, 'both') % no time or space constraint
        float_profs{fl} = all_prof(has_sensor & is_ocean & ...
            has_mode & has_press);
    else
        warning('no such setting for "outside": %s', outside)
        float_profs{fl} = [];
    end
 
    if isempty(float_profs{fl})
        float_ids(float_ids == good_float_ids(fl)) = [];
    end
end

float_profs(cellfun(@isempty, float_profs)) = [];

if min_num_prof
    has_num = cellfun(@length, float_profs) >= min_num_prof;
    float_ids = float_ids(has_num);
    float_profs = float_profs(has_num);
end
