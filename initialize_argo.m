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
%   can be used.
%
% AUTHORS: 
%   H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL),
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

global Settings Sprof Float;

addpath('auxil')
addpath(genpath('m_map'))
addpath(genpath('gsw'))

% The do_pause() function can be used in a driver script to
% halt the execution until the user presses ENTER.
% Set Settings.use_pause to 0 if you want to run everything without stopping.
use_desktop = desktop('-inuse');
Settings.use_pause = ~use_desktop;

% By default, actively running commands are described with output
% to the command windows. Set this to 0 to suppress this output.
% Values larger than 1 can be used to help in debugging.
Settings.verbose = 1;

% Maximum number of plots that can be created with one call to
% show_profiles etc.
% Increase this number if necessary, if you are sure that 
% your system can handle it
Settings.max_plots = 20;

% Profiles are stored in subdirectory 'Profiles'
Settings.prof_dir = './Profiles/';

% Index files are stored in subdirectory 'Index'
Settings.index_dir = './Index/'; 

% Create Index directory if needed
if ~check_dir(Settings.index_dir)
    error('Index directory could not be created')
end

% Create Profile directory if needed
if ~check_dir(Settings.prof_dir)
    error('Profile directory could not be created')
end

Settings.demo_float = 5904021;

% By default, don't update if files are less than 1 hour old
% alternative settings are 0 (don't update at all if files exist
% locally already) or 1 (always update)
Settings.update = 3600; % time is given in seconds

% default values for computation of mixed layer depth
Settings.temp_thresh = 0.2;
Settings.dens_thresh = 0.03;

% Default: try French GDAC before US GDAC
host_ifremer = 'https://data-argo.ifremer.fr/';
host_godae = 'https://usgodae.org/ftp/outgoing/argo/';
% Additional hosts could be added here
% Settings.hosts = {host_ifremer;host_godae}; % alternate order of hosts
Settings.hosts = {host_godae;host_ifremer};

Settings.avail_vars = {'PRES';'PSAL';'TEMP';'DOXY';'BBP';'BBP470';'BBP532';...
    'BBP700';'TURBIDITY';'CP';'CP660';'CHLA';'CDOM';'NITRATE';'BISULFIDE';...
    'PH_IN_SITU_TOTAL';'DOWN_IRRADIANCE';'DOWN_IRRADIANCE380';...
    'DOWN_IRRADIANCE412';'DOWN_IRRADIANCE443';'DOWN_IRRADIANCE490';...
    'DOWN_IRRADIANCE555';'DOWN_IRRADIANCE670';'UP_RADIANCE';...
    'UP_RADIANCE412';'UP_RADIANCE443';'UP_RADIANCE490';'UP_RADIANCE555';...
    'UP_RADIANCE';'UP_RADIANCE412';'UP_RADIANCE443';'UP_RADIANCE490';...
    'UP_RADIANCE555';'DOWNWELLING_PAR'};

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
sprof_date = H{2};
Sprof.date = datenum(sprof_date,'yyyymmddHHMMSS');
Sprof.lat  = H{3};
Sprof.lon  = H{4};
Sprof.ocean = H{5};
% column 5: ocean basin
% column 6: profiler type
% column 7: institution
Sprof.sens = H{8};
Sprof.data_mode = H{9};
Sprof.date_update = datenum(H{10}, 'yyyymmddHHMMSS');
Sprof.nprofs = length(H{1});

% Extract unique floats
Sprof.wmo = regexp(sprof_urls,'\d{7}','once','match');
[uwmo,ia] = unique(Sprof.wmo,'stable'); % keep list order
ulist = sprof_urls(ia);
dacs = regexp(ulist(:,1),'^\w+','once','match');
Sprof_fnames = regexprep(uwmo,'\d{7}','$0_Sprof.nc');
tmp = regexprep(ulist(:,1),'profiles.+','');
Sprof_fp = strcat(tmp,Sprof_fnames);

% Put per-float information into global struct Float
Float.file_path = Sprof_fp;
Float.file_name = Sprof_fnames;
Float.dac = dacs;
Float.wmoid = uwmo;
Float.nfloats = length(uwmo);
% range of profile indices per float
ia(end+1) = length(sprof_urls) + 1;
Float.prof_idx1 = ia(1:end-1);
Float.prof_idx2 = ia(2:end) - 1;
% use the update date of the last profile
Float.update = Sprof.date_update(Float.prof_idx2);

% Set up float/profile conversion matrix and profile-per-float IDs
%
% Details about the Sprof.fprofid array:
% It has N non-zero entries, where N is the total number of profiles that
% the Sprof index file contains, which corresponds to the number of lines
% in that file (minus 9, which is the number of its header lines).
% These entries are the overall indices of all profiles. 
% The values are the per-profile indices, i.e.,
% Sprof(i) = j  -> The i-th overall profile is the j-th profile for that
% particular float.
%
% The Sprof.p2f sparse matrix is used to select profiles of floats
% that have at least one profile that matches the given space and time
% constraints.
% Its dimensions are Sprof.nprofs x Float.nfloats.
% A vector with Sprof.nprofs entries that match the criteria (1=yes,0=no)
% multiplied by Sprof.p2f results in a vector that has positive values for
% all floats that have at least one matching profile, 0s for all floats
% that do not have any matching profiles. (The number is equal to the 
% number of profiles per float that match the given constraints.)
% Multiplying this result with the transpose of Sprof.p2f results in 
% a vector with Sprof.nprofs entries, positive for all profiles from all
% floats that have at least one matching profile, zeros for all others.
% See function select_profiles for the implementation of this 
% selection algorithm.
Sprof.p2f = sparse(Sprof.nprofs, Float.nfloats);
Sprof.fprofid = zeros(Sprof.nprofs, 1); % pre-allocate
for f = 1:Float.nfloats
    Sprof.p2f(Float.prof_idx1(f):Float.prof_idx2(f),f) = 1;
    Sprof.fprofid(Float.prof_idx1(f):Float.prof_idx2(f),1) = ...
        1:Float.prof_idx2(f) - Float.prof_idx1(f) + 1;
end

% Determine the availability of mapping functions
if ~isempty(which('geobasemap'))
    Settings.mapping = 'native';
elseif ~isempty(which('m_proj'))
    Settings.mapping = 'm_map';
else
    Settings.mapping = 'plain';
end
