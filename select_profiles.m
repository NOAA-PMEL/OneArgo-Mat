function [profiles,floats] = select_profiles(lon_lim,lat_lim,...
    start_date,end_date,varargin)
% select_profiles  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   [profiles,floats] = select_profiles(lon_lim,lat_lim,...
%                                       start_date,end_date,varargin)
%
% DESCRIPTION:
%   This function returns the indices of profiles and floats that match
%   the given criteria (spatial, temporal, sensor availability).
%   It does not download any files.
%
% INPUTS:
% lon_lim    : Longitude limits
% lat_lim    : Latitude limits
%            * Latitude and longitude limits can each be input as either a
%            two element vector (e.g., [LAT1 LAT2]) for maximum and minimum
%            limits or a four element vector (e.g., [LAT1 LAT2 LAT3 LAT4])
%            for quadrilateral vertices
%            * Longitude can be input in either the -180 to 180 degrees
%            format or 0 to 360 degrees format
%            * Both values can be '[]' to indicate the full range
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
%   profiles : array with the global indices of all matching profiles
%   floats   : array with the global indices of all matching floats
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

global Settings Sprof;

% set defaults
outside = 'none'; % if set, removes profiles outside time/space constraints
sensor = []; % default: use all available profiles
profiles = []; % will be assigned if successful
floats = []; % will be assigned if successful

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

% FORMAT LATITUDE AND LONGITUDE LIMITS
lon_lim(lon_lim>180) = lon_lim(lon_lim>180)-360; % Wrap lon to -180..180 deg
lon_lim(lon_lim<-180) = lon_lim(lon_lim<-180)+360;
% If only two elements are input, extend to obtain four vertices
if numel(lat_lim)==2
    latv = [lat_lim(1) lat_lim(1) lat_lim(2) lat_lim(2)];
else
    latv = lat_lim;
end
if numel(lon_lim)==2
    lonv = [lon_lim(1) lon_lim(2) lon_lim(2) lon_lim(1)];
else
    lonv = lon_lim;
end

% ADJUST INPUT DATES TO DATENUM FORMAT
dn1 = datenum(start_date); 
dn2 = datenum(end_date);

% GET INDEX OF PROFILES WITHIN USER-SPECIFIED GEOGRAPHIC POLYGON
if lonv(1) > lonv(2)
    lonv1 = [lonv(1) 180 180 lonv(1)];
    lonv2 = [-180 lonv(2) lonv(2) -180];
    inpoly = inpolygon(Sprof.lon,Sprof.lat,lonv1,latv) | ...
             inpolygon(Sprof.lon,Sprof.lat,lonv2,latv);
else
    inpoly = inpolygon(Sprof.lon,Sprof.lat,lonv,latv);
end

% Find index of dates that are within the time window
indate   = Sprof.date >= dn1 & Sprof.date <= dn2;

% SELECT BY SENSOR
if isempty(sensor)
    has_sensor = ones(size(indate)); % no sensor was selected
else
    has_sensor = contains(Sprof.sens, sensor);
end
all_prof = 1:length(indate);

if ~sum(has_sensor)
    warning('no data found for sensor %s', sensor);   
elseif strcmp(outside, 'none')
    profiles = all_prof(inpoly & indate & has_sensor);
else
    % identify all profiles of all those floats that have at least
    % one profile within given time and space constraints
    use_idx = ((inpoly' & indate') * Sprof.p2f) * Sprof.p2f';
    % now apply the given constraints
    if strcmp(outside, 'time') % must meet space constraint
        profiles = all_prof(inpoly & use_idx' & has_sensor);
    elseif strcmp(outside, 'space') % must meet time constraint
        profiles = all_prof(indate & use_idx' & has_sensor);
    elseif strcmp(outside, 'both') % no time or space constraint
        profiles = all_prof(use_idx' & has_sensor);
    else
        warning('no such setting for "outside": %s', outside)
    end
end

if ~isempty(profiles)
    floats = unique(Sprof.wmo(profiles));
end
