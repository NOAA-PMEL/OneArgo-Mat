function initialize_argo()
% initialize_argo  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   initialize_argo()
%
% DESCRIPTION:
%   This function defines standard settings and paths and downloads
%   index files. It must be called once before any other functions
%   can be used, either directly or indirectly by calling any of
%   the functions load_float_data, select_profiles, show_profiles,
%   show_sections, or show_trajectories.
%
% INPUT: None
%
% OUTPUT: None
%
% AUTHORS:
%   H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL), N. Buzby (UW)
%
% CITATION:
%   H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2022. OneArgo-Mat:
%   A MATLAB toolbox for accessing and visualizing Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.6588041
%
% LICENSE: oneargo_mat_license.m
%
% DATE: JUNE 1, 2022  (Version 1.0.1)

global Settings Prof Sprof Float;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGINNING OF SECTION WITH USER SPECIFIC OPTIONS
% this part of the function can be modified to meet specific needs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The do_pause() function can be used in a driver script to
% halt the execution until the user presses ENTER.
% Set Settings.use_pause to 0 if you want to run everything without stopping.
use_desktop = desktop('-inuse');
Settings.use_pause = ~use_desktop;

% By default, actively running commands are described with output
% to the command window. Set this to 0 to suppress this output.
% Values larger than 1 can be used to help in debugging.
Settings.verbose = 1;

% Maximum number of plots that can be created with one call to
% show_profiles etc.
% Increase this number if necessary, if you are sure that
% your system can handle it.
Settings.max_plots = 20;

% Profiles are stored in subdirectory 'Profiles'
Settings.prof_dir = './Profiles/';

% Index files are stored in subdirectory 'Index'
Settings.index_dir = './Index/';

% Meta files are stored in subdirectory 'Meta'
Settings.meta_dir = './Meta/';

% Tech files are stored in subdirectory 'Tech'
Settings.tech_dir = './Tech/';

% Traj files are stored in subdirectory 'Traj'
Settings.traj_dir = './Traj/';

Settings.demo_float = 5904021;

% By default, don't update index files if they are less than 1 hour old
% alternative settings are 0 (don't update at all if files exist
% locally already) or 1 (always update)
Settings.update = 3600; % time is given in seconds

% To ensure compatibility with the BGC-Argo-Mat toolbox
% (for calls to select_profiles that do not specify the
% type or BGC sensors), the default type can be set to
% 'bgc' here. 'phys' is the third available setting.
Settings.default_type = 'all';

% default values for computation of mixed layer depth
Settings.temp_thresh = 0.2;
Settings.dens_thresh = 0.03;

% default value for deviation from requested depth level for time series
% plots - if exceeded, a warning will be issued
Settings.depth_tol = 5;

% Settings.colormap = 'jet'; % uncomment and change as needed

% colors for profile plots ("range" is for individual profiles using the
% 'all' method and mean +- std.dev. in the 'mean' method)
Settings.color_var1_mean = [0, 0, 0]; % black
Settings.color_var1_range = [0.7 0.7 0.7]; % light gray
Settings.color_var2_mean = [0, 0, 1]; % blue
Settings.color_var2_range = [0.5, 0.75, 1]; % light blue

% color for estimated locations in trajectory plots
Settings.color_estim_loc = [0.7 0.7 0.7]; % light gray

% colors for data modes in trajectory plots:
% blue for R, yellow for A, green for D
Settings.traj_mode_colors = {[0, 0.4470, 0.7410]; ...
    [0.9290, 0.6940, 0.1250]; [0.4660, 0.6740, 0.1880]};

% amount of lon/lat padding in trajectory plots (in degrees)
Settings.pad_lon = 5;
Settings.pad_lat = 5;

% Default: try US GDAC before French GDAC
host_ifremer = 'https://data-argo.ifremer.fr/';
host_godae = 'https://usgodae.org/ftp/outgoing/argo/';
% Additional hosts could be added here
%Settings.hosts = {host_godae;host_ifremer};
% downloads from IFREMER are often faster than from GODAE
Settings.hosts = {host_ifremer;host_godae}; % alternate order of hosts

% Default: do not interpolate missing lon/lat values
Settings.interp_lonlat = 'no';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF SECTION WITH USER SPECIFIC OPTIONS
% the rest of this function should not be modified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add subdirectories with auxiliary functions to the path
filepath = fileparts(mfilename('fullpath'));
addpath([filepath, '/auxil'])
addpath(genpath([filepath, '/m_map']))
addpath(genpath([filepath, '/gsw']))

% Create subdirectories if needed
if ~check_dir(Settings.index_dir) || ~check_dir(Settings.prof_dir) || ...
        ~check_dir(Settings.meta_dir) || ~check_dir(Settings.tech_dir) || ...
        ~check_dir(Settings.traj_dir)
    error('Subdirectories could not be created')
end 

% Full set of available variables (but not all floats have all sensors)
% Additional sensors of existing types (e.g., DOXY2, BBP700_2) will
% be added below if they are found in the index file
Settings.avail_vars = {'PRES';'PSAL';'TEMP';'CNDC';'DOXY';'BBP';'BBP470';...
    'BBP532';'BBP700';'TURBIDITY';'CP';'CP660';'CHLA';'CDOM';'NITRATE';...
    'BISULFIDE';'PH_IN_SITU_TOTAL';'DOWN_IRRADIANCE';'DOWN_IRRADIANCE380';...
    'DOWN_IRRADIANCE412';'DOWN_IRRADIANCE443';'DOWN_IRRADIANCE490';...
    'DOWN_IRRADIANCE555';'DOWN_IRRADIANCE670';'DOWNWELLING_PAR';...
    'UP_RADIANCE';'UP_RADIANCE412';'UP_RADIANCE443';'UP_RADIANCE490';...
    'UP_RADIANCE555';};

% List of Data Assimilation Centers
Settings.dacs = {'aoml'; 'bodc'; 'coriolis'; 'csio'; 'csiro'; 'incois'; ...
    'jma'; 'kma'; 'kordi'; 'meds'};

% Download Sprof index file from GDAC to Index directory
sprof = 'argo_synthetic-profile_index.txt'; % file at GDAC
if ~download_index([sprof, '.gz'], 'Sprof')
    error('Sprof index file could not be downloaded')
end

% Extract information from Sprof index file
% NOTE that some quantities will be kept per float (struct Float):
% file_path, file_name, dac, params, wmoid, update
% Others will be kept per profile (struct Sprof):
% date, lat, lon, sens(ors), data_mode
initialize_sprof([Settings.index_dir, sprof]);

% Extract unique BGC floats
[uwmo_sprof,ia] = unique(Sprof.wmo,'stable'); % keep list order
Sprof_wmo = Sprof.wmo; % needed later if sprof_only is not empty
Sprof.wmo = str2double(Sprof.wmo);
ia(end+1) = length(Sprof.urls) + 1;
bgc_prof_idx1 = ia(1:end-1);
bgc_prof_idx2 = ia(2:end) - 1;

% Download prof index file from GDAC to Index directory
prof = 'ar_index_global_prof.txt'; % file at GDAC
if ~download_index([prof, '.gz'], 'prof')
    error('prof index file could not be downloaded')
end

% Extract information from prof index file
% NOTE that some quantities will be kept per float (struct Float):
% file_path, file_name, dac, params, wmoid, update
% Others will be kept per profile (struct Prof):
% date, lat, lon, sens(ors), data_mode
initialize_prof([Settings.index_dir, prof]);

% Extract unique floats
% note that older floats have 5-digit IDs
[uwmo_prof,ia2] = unique(Prof.wmo,'stable'); % keep list order
Prof_wmo = Prof.wmo; % needed later if sprof_only is not empty
Prof.wmo = str2double(Prof.wmo);
ulist = Prof.urls(ia2);
dacs = regexp(ulist(:),'^\w+','once','match');
prof_fnames = regexprep(uwmo_prof,'\d+','$0_prof.nc');
tmp = regexprep(ulist(:),'profiles.+','');
prof_fp = strcat(tmp,prof_fnames);

% are there any floats that have only Sprof files, but not prof files?
sprof_only = setdiff(uwmo_sprof, uwmo_prof);
if ~isempty(sprof_only)
    n_sprof_only = length(sprof_only);
    fprintf('The following %d floats have only Sprof files, not prof files:\n', ...
        n_sprof_only)
    disp(cell2mat(sprof_only))
    disp('') % empty line in command window
    idx_sprof_only = cellfun(@(x) find(strcmp(Sprof_wmo, x), 1), sprof_only);
    add_dacs = regexp(Sprof.urls(idx_sprof_only),'^\w+','once','match');
    % since the float type is 'bgc', prof will be replaced by Sprof later
    add_fnames = regexprep(sprof_only,'\d+','$0_prof.nc');
    tmp = regexprep(Sprof.urls(idx_sprof_only),'profiles.+','');
    add_fp = strcat(tmp, add_fnames);
    % add them to other variables so Float structure will have them
    prof_fp = [prof_fp; add_fp];
    prof_fnames = [prof_fnames; add_fnames];
    dacs = [dacs; add_dacs];
    wmo_sprof_only = str2double(sprof_only);
    idx_in_sprof = cell2mat(arrayfun(@(x) find(Sprof.wmo == x), ...
        wmo_sprof_only, 'UniformOutput', false));
    % add all fields from the "special" floats into the Prof structure
    Prof.urls = [Prof.urls; Sprof.urls(idx_in_sprof)];
    Prof.date = [Prof.date; Sprof.date(idx_in_sprof)];
    Prof.lat = [Prof.lat; Sprof.lat(idx_in_sprof)];
    Prof.lon = [Prof.lon; Sprof.lon(idx_in_sprof)];
    Prof.ocean = [Prof.ocean; Sprof.ocean(idx_in_sprof)];
    Prof.profiler = [Prof.profiler; Sprof.profiler(idx_in_sprof)];
    Prof.update = [Prof.update; Sprof.update(idx_in_sprof)];
    Prof.split_sens = [Prof.split_sens; Sprof.split_sens(idx_in_sprof)];
    Prof.wmo = [Prof.wmo; Sprof.wmo(idx_in_sprof)];
    Prof_wmo = [Prof_wmo; Sprof_wmo(idx_in_sprof)];
    [uwmo_prof,ia2] = unique(Prof_wmo,'stable');
end


% need to find out which floats are phys (in Prof only) and bgc (in Sprof)
is_uniq_bgc = ismember(uwmo_prof,uwmo_sprof);
nbgc = sum(is_uniq_bgc); % # of bgc floats (this may be revised later)

% determine index pointers from prof to Sprof files for all BGC floats
% (this needs to be done before the type is changed for those floats
% that are listed in Sprof index file but don't have BGC sensors)
bgc_idx_full = zeros(size(is_uniq_bgc));
bgc_idx_full(is_uniq_bgc) = 1:nbgc;

% Put per-float information into global struct Float
Float.file_path = prof_fp;
Float.file_name = prof_fnames;
Float.dac = dacs;
Float.wmoid = str2double(uwmo_prof);
Float.nfloats = length(uwmo_prof);
% range of profile indices per float, referring to the Prof struct
ia2(end+1) = length(Prof.urls) + 1;
Float.prof_idx1 = ia2(1:end-1);
Float.prof_idx2 = ia2(2:end) - 1;
Float.profiler = Prof.profiler(Float.prof_idx1);
Float.type = cell(Float.nfloats, 1);
Float.type(cellfun(@isempty, Float.type)) = {'phys'};
Float.type(is_uniq_bgc) = {'bgc'};

% determine types of sensors/variables that are present for some and for
% all profiles of any given float; also re-flag floats from BGC to phys
% if they don't have any BGC variables available
Float.min_sens = cell(Float.nfloats, 1); % pre-allocate cell arrays
Float.max_sens = cell(Float.nfloats, 1);
len_sens = cellfun(@length, Sprof.sens);
count = 0;
index_bgc = 0;
is_true_bgc = ones(length(bgc_prof_idx1), 1);
count_bgc = 0; % used for index_full_to_bgc
index_full_to_bgc = nan(Float.nfloats, 1); % used for finding update dates
for f = 1:Float.nfloats
    if strcmp(Float.type{f}, 'phys')
        % ar_index_global_prof.txt does not contain information about
        % the available sensors per profile
        if Float.profiler(f) == 845
            Float.min_sens{f} = {'PRES';'TEMP'};
        else
            Float.min_sens{f} = {'PRES';'TEMP';'PSAL'};
        end
        Float.max_sens{f} = Float.min_sens{f};
    else % BGC float
        index_bgc = index_bgc + 1;
        f2 = bgc_idx_full(f);
        [~, idx1] = min(len_sens(bgc_prof_idx1(f2):bgc_prof_idx2(f2)));
        [~, idx2] = max(len_sens(bgc_prof_idx1(f2):bgc_prof_idx2(f2)));
        % assumption: the shortest string has sensors that are shared among all
        % profiles and the longest string has the union of all available sensors
        Float.min_sens{f} = Sprof.split_sens{bgc_prof_idx1(f2) + idx1 - 1};
        Float.max_sens{f} = Sprof.split_sens{bgc_prof_idx1(f2) + idx2 - 1};
        % if there are no profiles for this float with any BGC sensors
        % set its type to 'phys'
        bgc_sensors = Float.max_sens{f};
        bgc_sensors(ismember(bgc_sensors, ...
            {'PRES';'TEMP';'PSAL';'CNDC'})) = [];
        if isempty(bgc_sensors)
            Float.type{f} = 'phys';
            count = count + 1;
            is_true_bgc(index_bgc) = 0;
        else
            count_bgc = count_bgc + 1;
            index_full_to_bgc(f) = count_bgc;
        end
    end
end
fprintf('Note: %d floats from Sprof index file do not have BGC sensors\n', ...
    count);
% assign index of first and last profile to true BGC floats only, referring
% to the indices within the Sprof struct
% these variables should never be used for non-BGC floats, so their value
% of 0 serves as a flag that would result in an out-of-bounds error
Float.bgc_prof_idx1 = zeros(size(Float.prof_idx1));
Float.bgc_prof_idx2 = Float.bgc_prof_idx1;
bgc_prof_idx1 = bgc_prof_idx1(is_true_bgc == 1);
bgc_prof_idx2 = bgc_prof_idx2(is_true_bgc == 1);
Float.bgc_prof_idx1(strcmp(Float.type, 'bgc')) = bgc_prof_idx1;
Float.bgc_prof_idx2(strcmp(Float.type, 'bgc')) = bgc_prof_idx2;
% add this check in case that files without prof files do not have
% any BGC sensors, so that they will be classified as 'phys' floats
% but Sprof files need to be read anyway
Float.has_prof_file = ones(length(Float.type), 1);
if ~isempty(sprof_only)
    idx_in_float = cell2mat(arrayfun(@(x) find(Float.wmoid == x), ...
        wmo_sprof_only, 'UniformOutput', false));
    Float.has_prof_file(idx_in_float) = 0;
end

% for all "true" BGC floats, i.e., those that are listed in the Sprof
% index file and have more than pTS sensors, Sprof files will be used
% instead of prof files
idx_bgc = strcmp(Float.type, 'bgc');
fprintf('%d "true" BGC floats were found\n', sum(idx_bgc));
idx_phys = strcmp(Float.type, 'phys');
fprintf('%d core and deep floats were found\n', sum(idx_phys));
Float.file_path(strcmp(Float.type, 'bgc') | ~Float.has_prof_file) = ...
    cellfun(@(x) strrep(x, 'prof', 'Sprof'), ...
    Float.file_path(strcmp(Float.type, 'bgc') | ~Float.has_prof_file), ...
    'UniformOutput', false);
Float.file_name(strcmp(Float.type, 'bgc') | ~Float.has_prof_file) = ...
    cellfun(@(x) strrep(x, 'prof', 'Sprof'), ...
    Float.file_name(strcmp(Float.type, 'bgc') | ~Float.has_prof_file), ...
    'UniformOutput', false);

% use the most recent update date across profiles for any given float
for f = 1:length(Float.prof_idx1)
    if strcmp(Float.type{f}, 'bgc')
        flt = index_full_to_bgc(f);
        update_dates = str2double(Sprof.update(bgc_prof_idx1(flt):bgc_prof_idx2(flt)));
        [~,idx] = max(update_dates);
        Float.update{f} = Sprof.update{bgc_prof_idx1(flt) + idx - 1};
    else
        update_dates = str2double(Prof.update(Float.prof_idx1(f):Float.prof_idx2(f)));
        [~,idx] = max(update_dates);
        Float.update{f} = Prof.update{Float.prof_idx1(f) + idx - 1};
    end
end

% Download meta index file from GDAC to Index directory
% Since it is rather small, download the uncompressed file directly
meta = 'ar_index_global_meta.txt';
if ~download_index(meta, 'meta')
    error('meta index file could not be downloaded')
end

% Extract information from meta index file
initialize_meta([Settings.index_dir, meta]);

% Download tech index file from GDAC to Index directory
% Since it is rather small, download the uncompressed file directly
tech = 'ar_index_global_tech.txt';
if ~download_index(tech, 'tech')
    error('tech index file could not be downloaded')
end

% Extract information from tech index file
initialize_tech([Settings.index_dir, tech]);

% Download traj index file from GDAC to Index directory
% Since it is rather small, download the uncompressed file directly
traj = 'ar_index_global_traj.txt';
if ~download_index(traj, 'traj')
    error('traj index file could not be downloaded')
end

% Extract information from traj index file
initialize_traj([Settings.index_dir, traj]);

% Determine the availability of mapping functions
if ~isempty(which('geobasemap'))
    Settings.mapping = 'native';
elseif ~isempty(which('m_proj'))
    Settings.mapping = 'm_map';
else
    Settings.mapping = 'plain';
end
