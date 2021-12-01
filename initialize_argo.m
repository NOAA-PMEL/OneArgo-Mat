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
Settings.avail_vars = {'PRES';'PSAL';'TEMP';'DOXY';'BBP';'BBP470';'BBP532';...
    'BBP700';'TURBIDITY';'CP';'CP660';'CHLA';'CDOM';'NITRATE';'BISULFIDE';...
    'PH_IN_SITU_TOTAL';'DOWN_IRRADIANCE';'DOWN_IRRADIANCE380';...
    'DOWN_IRRADIANCE412';'DOWN_IRRADIANCE443';'DOWN_IRRADIANCE490';...
    'DOWN_IRRADIANCE555';'DOWN_IRRADIANCE670';'UP_RADIANCE';...
    'UP_RADIANCE412';'UP_RADIANCE443';'UP_RADIANCE490';'UP_RADIANCE555';...
    'UP_RADIANCE';'UP_RADIANCE412';'UP_RADIANCE443';'UP_RADIANCE490';...
    'UP_RADIANCE555';'DOWNWELLING_PAR';'DOXY2';'DOXY3'};

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
Sprof.data_mode = H{9};
Sprof.date_update = H{10};

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

% Determine the availability of mapping functions
if ~isempty(which('geobasemap'))
    Settings.mapping = 'native';
elseif ~isempty(which('m_proj'))
    Settings.mapping = 'm_map';
else
    Settings.mapping = 'plain';
end
