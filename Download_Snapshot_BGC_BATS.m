%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Download_Snapshot_BGC_BATS.m
% Driver routine for the MATLAB toolbox for accessing OneArgo float data.
%
% Initially written for the FAIR workshop BGC-Argo Tutorial
% April 21-25, 2025
%
% Demonstrates the downloading of a snapshot of BGC-Argo data for the
% vicinity of the BATS location. Temporarily, over 10 GB of local disk
% space is needed. After unpacking and deleting unneeded files, only
% 300 MB of local disk space is needed for the snapshot files, if the
% first option is chosen. If all Sprof files from the snapshot are kept,
% 14 GB of disk space are needed for the snapshot directory.
%
% TUTORIAL AUTHOR:
%   H. Frenzel (UW-CICOES)
%
% OneArgo-Mat AUTHORS:
%   J. Sharp and H. Frenzel (UW-CICOES), A. Fassbender (NOAA-PMEL), N. Buzby (UW)
%
% CITATION:
%   H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2025. OneArgo-Mat:
%   A MATLAB toolbox for accessing and visualizing Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.6588041
%
% LICENSE: oneargo_mat_license.m
%
% DATE: APRIL 16, 2025  (Version 1.1.0)

global Settings;

bats_lon = -64.17;
bats_lat = 31.67;
pad_xy = 5; % find floats within 5 degrees in longitude and latitude

% find BGC floats from any time
bats_floats = select_profiles([bats_lon-pad_xy, bats_lon+pad_xy], ...
    [bats_lat-pad_xy, bats_lat+pad_xy], [], [], 'type', 'bgc');

fprintf('\n%d BGC floats within %d degrees of the BATS site were found.\n', ...
    length(bats_floats), pad_xy)

disp('Now the April 2025 snapshot of BGC Argo data will be downloaded.')
disp('This will take at least 10 minutes.')
while 1
    disp('You have two options:')
    disp('1: Only the Sprof files for floats that were found near BATS will be kept.')
    disp('2: All Sprof files will be kept.')
    disp('Select this option if you want to run the full tutorial with snapshot data.')
    result = input('Your choice: ', 's');
    if strcmp(result, '1') || strcmp(result, '2')
        break
    else
        disp('Please enter "1" or "2"!')
    end
end
choice = uint8(str2double(result));
if choice == 1
    vargin = {'floats', bats_floats};
else
    vargin = {};
end

download_snapshot('type', 'bgc', 'date', 202504, vargin{:})

Settings.default_type = 'bgc';
Settings.use_snapshots = 202504;

% use all files from the snapshot, not the GDAC
Settings.prof_dir = [Settings.snap_dir, Settings.snap_path, 'Profiles/'];
Settings.index_dir = [Settings.snap_dir, Settings.snap_path, 'Index/'];
Settings.meta_dir = [Settings.snap_dir, Settings.snap_path, 'Meta/'];
Settings.tech_dir = [Settings.snap_dir, Settings.snap_path, 'Tech/'];
Settings.traj_dir = [Settings.snap_dir, Settings.snap_path, 'Traj/'];

% there are two floats that are listed in argo_synthetic-profile_index.txt,
% but do not have any BGC data
bats_floats(bats_floats == 4902441) = [];
bats_floats(bats_floats == 4902524) = [];

% one sample plot using the snapshot files
show_trajectories(bats_floats, 'color', 'multiple', 'legend', 'n');

%%
fprintf('\n')
disp('To keep working with these snapshot data, make the following changes')
disp('in initialize_argo.m:')
disp('Line 94: set Settings.default_type to ''bgc''.')
disp('Comment out line 108 ("Settings.use_snapshots = 0;").')
disp('Uncomment line 112 ("Settings.use_snapshots = 202504;").')

if choice == 2
    fprintf('\n')
    disp('Afterwards, you can run FAIR_Tutorial again, using the snapshot data.');
end

