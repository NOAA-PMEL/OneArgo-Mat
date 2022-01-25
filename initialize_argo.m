function initialize_argo()
% initialize_argo  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
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
% AUTHORS:
%   H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL), N. Buzby (UW),
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

global Settings Sprof Float;

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

Settings.demo_float = 5904021;

% By default, don't update if files are less than 1 hour old
% alternative settings are 0 (don't update at all if files exist
% locally already) or 1 (always update)
Settings.update = 3600; % time is given in seconds

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

% amount of lon/lat padding in trajectory plots (in degrees)
Settings.pad_lon = 5;
Settings.pad_lat = 5;

% Default: try US GDAC before French GDAC
host_ifremer = 'https://data-argo.ifremer.fr/';
host_godae = 'https://usgodae.org/ftp/outgoing/argo/';
% Additional hosts could be added here
Settings.hosts = {host_godae;host_ifremer};
% Settings.hosts = {host_ifremer;host_godae}; % alternate order of hosts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF SECTION WITH USER SPECIFIC OPTIONS
% the rest of this function should not be modified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add subdirectories with auxiliary functions to the path
filepath = fileparts(mfilename('fullpath'));
addpath([filepath, '/auxil'])
addpath(genpath([filepath, '/m_map']))
addpath(genpath([filepath, '/gsw']))

% Create Index directory if needed
if ~check_dir(Settings.index_dir)
    error('Index directory could not be created')
end

% Create Profile directory if needed
if ~check_dir(Settings.prof_dir)
    error('Profile directory could not be created')
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

% Write Sprof index file from GDAC to Index directory
sprof = 'argo_synthetic-profile_index.txt'; % file used locally
sprof_gz = [sprof, '.gz']; % file that is available at GDAC
Settings.dest_path_sprof = [Settings.index_dir, sprof];
dest_path_sprof_gz = [Settings.index_dir, sprof_gz];
if do_download(dest_path_sprof_gz)
    if Settings.verbose
        disp('Sprof index file will now be downloaded.')
        disp('Depending on your internet connection, this may take a while.')
    end
    if ~try_download(sprof_gz, dest_path_sprof_gz)
        error('Sprof index file could not be downloaded')
    end
    gunzip(dest_path_sprof_gz)
elseif ~exist(Settings.dest_path_sprof, 'file')
    gunzip(dest_path_sprof_gz)
end

% Extract information from Sprof index file
% NOTE that some quantities will be kept per float (struct Float):
% file_path, file_name, dac, params, wmoid, update
% Others will be kept per profile (struct Sprof):
% date, lat, lon, sens(ors), data_mode
fid = fopen(Settings.dest_path_sprof);
H = textscan(fid,'%s %s %f %f %s %d %s %s %s %s','headerlines',9,...
    'delimiter',',','whitespace','');
fclose(fid);
sprof_urls = H{1};
Sprof.date = H{2};
Sprof.lat  = H{3};
Sprof.lon  = H{4};
Sprof.ocean = H{5};
% column 6: profiler type
% column 7: institution
Sprof.sens = H{8};
Sprof.split_sens = cellfun(@split, Sprof.sens, 'UniformOutput', false);
Sprof.data_mode = H{9};
Sprof.date_update = H{10};

% check for additional (e.g., DOXY2) and unknown sensors
for i = 1:length(Sprof.sens)
    sensors = split(Sprof.sens{i});
    if ~all(ismember(sensors, Settings.avail_vars))
        unknown_sensors = sensors(~ismember(sensors, Settings.avail_vars));
        for s = 1:length(unknown_sensors)
            main_sensor = get_sensor_number(unknown_sensors{s});
            if isempty(main_sensor)
                warning('unknown sensor in index file: %s', unknown_sensors{s})
            else
                Settings.avail_vars{end+1} = unknown_sensors{s};
            end
        end
    end
end

% Extract unique floats
Sprof.wmo = regexp(sprof_urls,'\d{7}','once','match');
[uwmo,ia] = unique(Sprof.wmo,'stable'); % keep list order
Sprof.wmo = str2double(Sprof.wmo);
ulist = sprof_urls(ia);
dacs = regexp(ulist(:,1),'^\w+','once','match');
Sprof_fnames = regexprep(uwmo,'\d{7}','$0_Sprof.nc');
tmp = regexprep(ulist(:,1),'profiles.+','');
Sprof_fp = strcat(tmp,Sprof_fnames);

% Put per-float information into global struct Float
Float.file_path = Sprof_fp;
Float.file_name = Sprof_fnames;
Float.dac = dacs;
Float.wmoid = str2double(uwmo);
Float.nfloats = length(uwmo);
% range of profile indices per float
ia(end+1) = length(sprof_urls) + 1;
Float.prof_idx1 = ia(1:end-1);
Float.prof_idx2 = ia(2:end) - 1;
% use the update date of the last profile
Float.update = Sprof.date_update(Float.prof_idx2);

% determine sensor availability by float
nfloats = length(Float.wmoid);
Float.min_sens = cell(nfloats, 1);
Float.max_sens = cell(nfloats, 1);
len_sens = cellfun(@length, Sprof.sens);
for f = 1:nfloats
    [~, idx1] = min(len_sens(Float.prof_idx1(f):Float.prof_idx2(f)));
    [~, idx2] = max(len_sens(Float.prof_idx1(f):Float.prof_idx2(f)));
    % assumption: the shortest string has sensors that are shared among all
    % profiles and the longest string has the union of all available sensors
    Float.min_sens{f} = Sprof.split_sens{Float.prof_idx1(f)+idx1-1};
    Float.max_sens{f} = Sprof.split_sens{Float.prof_idx1(f)+idx2-1};
end

% Determine the availability of mapping functions
if ~isempty(which('geobasemap'))
    Settings.mapping = 'native';
elseif ~isempty(which('m_proj'))
    Settings.mapping = 'm_map';
else
    Settings.mapping = 'plain';
end
