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
%           sensor type. Available are: DOXY, CHLA, BBP700, 
%           PH_IN_SITU_TOTAL, NITRATE, DOWN_IRRADIANCE380,
%           DOWN_IRRADIANCE412, DOWN_IRRADIANCE490, DOWNWELLING_PAR
%           (Currently, only one sensor type can be selected.)
%
% OUTPUTS:
%   float_ids   : array with the WMO IDs of all matching floats
%   float_profs : cell array with the per-float indices of all matching profiles
%
% AUTHORS: 
%   J. Sharp, H. Frenzel, A. Fassbender (NOAA-PMEL),
%   J. Plant, T. Maurer, Y. Takeshita (MBARI), D. Nicholson (WHOI),
%   and A. Gray (UW)
%
% CITATION:
%   H. Frenzel*, J. Sharp*, A. Fassbender, J. Plant, T. Maurer,
%   Y. Takeshita, D. Nicholson, A. Gray, 2021. BGC-Argo-Mat: A MATLAB
%   toolbox for accessing and visualizing Biogeochemical Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.4971318.
%   (*These authors contributed equally to the code.)
%
% LICENSE: bgc_argo_mat_license.m
%
% DATE: June 15, 2021

global Float Settings Sprof;

% set defaults
outside = 'none'; % if set, removes profiles outside time/space constraints
sensor = []; % default: use all available profiles

% parse optional arguments
for i = 1:2:length(varargin)
    if strcmpi(varargin{i}, 'outside')
        outside = varargin{i+1};
    elseif strcmpi(varargin{i}, 'sensor')
        sensor = varargin{i+1};
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
if isempty(inpoly)
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
if ~sum(has_sensor)
    warning('no data found for sensor %s', sensor);
end

all_prof = 1:length(indate);
profiles = all_prof(inpoly & indate & has_sensor);
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
    n_prof = get_dims(filename);
    fl_idx = find(str2num(cell2mat(Float.wmoid)) == good_float_ids(fl));
    n_prof_exp = Float.prof_idx2(fl_idx) - Float.prof_idx1(fl_idx) + 1;
    if n_prof_exp > n_prof
        warning(['The index file lists %d profiles for float %d, ', ...
            'but the Sprof file has only %d profiles.'], ...
            n_prof_exp, good_float_ids(fl), n_prof)
    end
    lon = ncread(filename, 'LONGITUDE');
    lat = ncread(filename, 'LATITUDE');
    juld = ncread(filename, 'JULD');
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
 
    % now apply the given constraints
    all_prof = 1:length(inpoly);
    if strcmp(outside, 'none')
        float_profs{fl} = all_prof(inpoly & indate & has_sensor);
    elseif strcmp(outside, 'time') % must meet space constraint
        float_profs{fl} = all_prof(inpoly & has_sensor);
    elseif strcmp(outside, 'space') % must meet time constraint
        float_profs{fl} = all_prof(indate & has_sensor);
    elseif strcmp(outside, 'both') % no time or space constraint
        float_profs{fl} = all_prof(has_sensor);
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
