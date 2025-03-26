function success = check_existing_snapshot()
% check_existing_snapshot  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   success = check_existing_snapshot(f)
%
% DESCRIPTION:
%   This function checks if files for the specified snapshot exist
%   already. It is called by determine_snapshot if the SEANOE website
%   cannot be reached.
%   Note that only the existence of a few directories and files
%   is checked.
%
% INPUTS:
%   None
%
% OUTPUT:
%   success   : 1 = snapshot files exist; 0 = snapshot fles are missing
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

global Settings;

success = 0; % change it once all tests have passed

if ~Settings.use_snapshots
    success = 1; % a snapshot is not requested
    return
end

if Settings.use_snapshots == 1
    % If the most recent snapshot is requested, we would normally check
    % seanoe.org what the most recent available snapshot is.
    % Since the website is presumably down, we will make an assumption
    % based on the current date and alert the user appropriately.
    snap_month = month(today);
    snap_year = year(today);
    if day(today) < 15
        % assume that snapshots are made available on the 15th for the
        % current month; before that, use the previous month
        snap_month = snap_month - 1;
        if snap_month == 0
            snap_month = 12;
            snap_year = snap_year - 1;
        end
    end
    curr_yr_mth = snap_year * 100 + snap_month;
    warning(['You requested the latest available snapshot, but ', ...
        'seanoe.org cannot be reached at this point.'])
    fprintf('Based on today''s date (%s), we assume that the most recent\n', ...
        datestr(today))
    fprintf('available snapshot is from %d.\n', curr_yr_mth)
    yr_mth = 0;
    while ~yr_mth
        fprintf('Please enter "1" if %d is the correct year and month ', ...
            curr_yr_mth)
        yr_mth = input('or enter the correct year and month in YYYYMM format: ');
        if yr_mth ~= 1 && (yr_mth < 201500 || yr_mth > curr_yr_mth || ...
                mod(yr_mth, 100) == 0 || mod(yr_mth, 100) > 12)
            disp('Please enter "1" or the desired year and month in the correct format!')
            yr_mth = 0;
        end
    end
    if yr_mth == 1
        snap_year_month = curr_yr_mth;
    else
        snap_year_month = yr_mth;
    end
else % a specific snapshot was requested
    snap_year_month = Settings.use_snapshots;
end
if strcmp(Settings.default_type, 'bgc')
    this_snap_dir = sprintf('%d-BgcArgoSprof', snap_year_month);
else % Settings.default_type is 'phys' or 'all'
    this_snap_dir = sprintf('%d-ArgoData', snap_year_month);
end
snap_path = [Settings.snap_dir, this_snap_dir];

dirs_to_check = {snap_path; [snap_path, '/Index']; [snap_path, '/Profiles']};
for d = 1:length(dirs_to_check)
    if exist(dirs_to_check{d}, 'dir') ~= 7
        warning('Directory "%s" not found!', dirs_to_check{d})
        return
    end
end

files_to_check = {[snap_path, '/Index/argo_synthetic-profile_index.txt']};
if ~strcmp(Settings.default_type, 'bgc')
    files_to_check{end+1} = [snap_path, '/Index/ar_index_global_prof.txt'];
end
for f = 1:length(files_to_check)
    if exist(files_to_check{f}, 'file') ~= 2
        warning('File "%s" not found!', files_to_check{f})
        return
    end
end

prof_files = dir([snap_path, '/Profiles/*prof.nc']);
if isempty(prof_files)
    warning('No profile files found in %s/Profiles!', snap_path)
    return
end

Settings.snap_path = snap_path;

success = 1; % all checks passed
