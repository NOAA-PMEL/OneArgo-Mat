function download_snapshot(varargin)
% download_snapshot  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   download_snapshot(varargin)
%
% DESCRIPTION:
%   This function downloads and unpacks one snapshot of Argo data.
%   By default, Settings.default_type and Settings.use_snapshot
%   values are considered to determine which snapshot will be downloaded.
%   Both settings can be overridden with optional inputs type and snap_date.
%
% INPUT: None
%
% OPTIONAL INPUTS  (key,value pairs):
%   'type', type: 'all' for all Argo data, 'phys' for core/deep only,
%         or 'bgc' for BGC Sprof files only
%   'dac, dac: keep only files from the specified DAC(s), e.g., 'AOML'
%         or {'jma';'coriolis'}
%   'date', date: 1 denotes the most recent available snapshot,
%         YYYYMM describes a particular snapshot (e.g., 202309 for
%         September 2023) - specified as integer, not as string;
%         if 0, nothing will be downloaded
%   'clobber', clobber: if 0 (default), do not download a snapshot file if
%         at least some of its contents exist locally already;
%         if 1, download snapshot in any case
%   'floats', floats: keep only files for the float(s) with the specified
%         WMO ID(s)
%   'keep', keep: if 0 (default), all files that are not used by this
%         toolbox are deleted; if 1, tech, meta, and traj files are kept
%         from full snapshots;
%         if 2, all files are kept (full tarball, intermediate tarballs
%         and individual profile netcdf files for full snapshots,
%         aux and geo subdirectories) -
%         note that this requires a lot of disk space, esp. for full snapshots
%         (If 'dac' or 'floats' options are used, non-matching float
%         profiles will be deleted even if keep is set to 2.)
%   'verbose', verbose: if 1, show more information about what is being
%         done; default: 0
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

global Settings;

% make sure Settings is initialized
if isempty(Settings)
    initialize_argo();
end

% set default values
snap_type = Settings.default_type;
snap_date = Settings.use_snapshots;
dac = [];
floats = [];
clobber = 0;
keep = 0;
verbose = 0;

% parse optional arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'type')
        snap_type = varargin{i+1};
    elseif strcmpi(varargin{i}, 'date')
        snap_date = varargin{i+1};
    elseif strcmpi(varargin{i}, 'dac')
        dac = lower(varargin{i+1});
    elseif strcmpi(varargin{i}, 'floats')
        floats = lower(varargin{i+1});
    elseif strcmpi(varargin{i}, 'clobber')
        clobber = varargin{i+1};
    elseif strcmpi(varargin{i}, 'keep')
        keep = varargin{i+1};
    elseif strcmpi(varargin{i}, 'verbose')
        verbose = varargin{i+1};
    end
end

if ~snap_date
    if ~isempty(varargin)
        % no need to issue a message if it was called from initialize_argo
        % without arguments
        disp('Settings.use_snapshots is set to 0, no snapshot will be downloaded.')
    end
    return
end
if snap_date > 1 && (snap_date < 201212 || snap_date > 203800)
    warning('download_snapshot: Wrong format of snap_date, must be 1 or YYYYMM')
    return
end

% set search patterns for parsing of web page
if strcmp(snap_type, 'all') || strcmp(snap_type, 'phys')
    pattern = 'Global GDAC Argo data files (';
    pat_month = 'Global GDAC Argo data files \(20\d{2}-\d{2}-\d{2} snapshot';
elseif strcmp(snap_type, 'bgc')
    pattern = 'BGC Sprof data files (';
    pat_month = 'BGC Sprof data files \(20\d{2}-\d{2}-\d{2} snapshot\)';
else
    fprintf('No such type of snapshot file: "%s"\n', snap_type);
    return
end

% always use cells for DAC(s); convert char array to cell
if ischar(dac)
    dac = cellstr(dac);
end

if isempty(dac) && ~isempty(floats)
    % select only the DACs that handle the specified floats
    dac = list_dacs(floats);
end

mth_idx = length(pattern) + 1; % first character of YYYY-MM

% download the web page that contains the links to all snapshots
page = webread('https://www.seanoe.org/data/00311/42182/');
idx = strfind(page, pattern);
if isempty(idx)
    fprintf('The snapshot web page could not be read properly.\n')
    return
end

% determine which snapshots are available
snap_month = [];
snap_size = [];
snap_url = {};
count = 0;
for i = 1:length(idx)
    this_link = page(idx(i):idx(i)+400);
    match_month = regexp(this_link, pat_month, 'match', 'once');
    pat_url = '"fileUrl":"http[\w/:.]+"';
    match_url = regexp(this_link, pat_url, 'match', 'once');
    pat_size = '"size":\d+';
    match_size = regexp(this_link, pat_size, 'match', 'once');
    % note that some snapshots are only available on demand, for those,
    % there is no fileUrl entry
    if ~isempty(match_month) && ~isempty(match_url)
        count = count + 1;
        year = str2double(match_month(mth_idx:mth_idx+3));
        month = str2double(match_month(mth_idx+5:mth_idx+6));
        snap_month(count) = 100 * year + month;
        snap_url{count} = match_url(12:end-1);
        if ~isempty(match_size)
            snap_size(count) = uint64(str2double(match_size(8:end)));
        else
            snap_size(count) = -1;
        end
    end
end

if ~count
    fprintf('No matching snapshots were found.\n')
    return
end

if snap_date == 1
    % use the most recent snapshot, need to sort them by snap_date
    [~,isort] = sort(snap_month,'descend'); % most recent will be first
    isnap = isort(1);
    snap_date = snap_month(isnap); % assignment makes it easier going forward
else
    % check if the specified snapshot exists
    [mini,isnap] = min(abs(snap_month - snap_date));
    if mini > 0
        fprintf('No snapshot found for %d\n', snap_date)
        fprintf('The nearest available snapshot is %d\n', snap_month(isnap))
        return
    end
end

if strcmp(snap_type, 'all') || strcmp(snap_type, 'phys')
    snap_path = sprintf('%s%d-ArgoData/', Settings.snap_dir, snap_date);
else
    snap_path = sprintf('%s%d-BgcArgoSprof/', Settings.snap_dir, snap_date);
end
Settings.snap_path = snap_path;

if exist(snap_path, 'dir') && ~clobber
    if ~isempty(varargin)
        fprintf('Snapshot for %d was downloaded already\n', snap_date);
    end
    return
end

% this is the most commonly used format for file names:
filename = regexp(snap_url{isnap}, '\d+\.tar.gz', 'match', 'once');
if isempty(filename) % alternate file name format
    filename = regexp(snap_url{isnap}, '\d+\.tgz', 'match', 'once');
end
if isempty(filename) % this should not happen
    fprintf('unexpected file name for snapshot: "%s"', snap_url{isnap})
    return
end

fprintf('Starting download of "%s" now\n', filename)
fprintf('Its file size is %.1f GB - this may take a while.\n', ...
    snap_size(isnap)*1e-9)
websave(filename, snap_url{isnap});
if ~exist(filename, 'file')
    fprintf('"%s" could not be downloaded\n', filename)
    return
end
% check if size of downloaded file matches the expected size
file = dir(filename);
file_size = file.bytes;
if file_size ~= snap_size(isnap)
    fprintf('"%s" was not downloaded correctly.\n', filename)
    fprintf('Expected file size: %d bytes\n', snap_size(isnap))
    fprintf('Actual file size:   %d bytes\n', file_size)
    return
end

try
    fprintf('Untarring the snapshot, this may take a few minutes...')
    untar(filename, Settings.snap_dir);
    fprintf('done!\n')
catch
    fprintf('error! Aborting...\n')
    return
end

% the original file is still present after untarring
if keep < 2
    delete(filename); % delete it to free up disk space
end

% Individual (S)prof files are sorted by dac subdirectory, which is
% a different approach from this toolbox's normal way of organizing files.
% Therefore, files will be reorganized.
snap_path_dac = [snap_path, 'dac'];

% In the case of full snapshots, there are core and bgc tarballs for
% each DAC, which must be unpacked first.
if ~strcmp(snap_type, 'bgc')
    % contents of aux and geo directory trees are not used by this toolbox
    if keep < 2
        rmdir([snap_path, '/aux'], 's');
        rmdir([snap_path, '/geo'], 's');
    end
    tarballs = dir([snap_path_dac, '/*.tar.gz']);
    fn_tar = {tarballs.name};
    for i = 1:length(fn_tar)
        if strcmp(snap_type, 'phys') && ...
                endsWith(fn_tar{i}, 'bgc.tar.gz') && keep < 2
            if verbose
                fprintf('Deleting "%s"\n', fn_tar{i})
            end
            delete([snap_path_dac, '/', fn_tar{i}])
            continue
        end
        if ~isempty(dac) && keep < 2
            keep_tar = 0;
            for d = 1:length(dac)
                if startsWith(fn_tar{i}, dac{d})
                    keep_tar = 1;
                    break;
                end
            end
            if ~keep_tar
                if verbose
                    fprintf('Deleting "%s"\n', fn_tar{i})
                end
                delete([snap_path_dac, '/', fn_tar{i}])
                continue
            end
        end
        file = dir([snap_path_dac, '/', fn_tar{i}]);
        file_size = file.bytes;
        if file_size > 1e9
            fprintf('Untarring "%s" (%.0f GB) may take a few minutes...\n', ...
                fn_tar{i}, 1e-9*file_size);
        elseif verbose
            fprintf('Untarring "%s"...\n', fn_tar{i})
        end
        try
            untar([snap_path_dac, '/', fn_tar{i}], snap_path_dac);
            fprintf('done!\n')
            if keep < 2
                delete([snap_path_dac, '/', fn_tar{i}])
            end
        catch
            fprintf('Error during untar!\n')
        end
        % the full snapshot archives include individual profile files,
        % which are not used by this toolbox
        if keep < 2
            match_dac = regexp(fn_tar{i}, '^[a-z]+_', 'match', 'once');
            this_dac = match_dac(1:end-1); % without the trailing _
            dac_contents = dir([snap_path_dac, '/', this_dac]);
            for j = 1:length(dac_contents)
                if dac_contents(j).isdir && ...
                        ~strcmp(dac_contents(j).name, '.') && ...
                        ~strcmp(dac_contents(j).name, '..')
                    dir_prof = sprintf('%s/%s/%s/profiles', snap_path_dac, ...
                        this_dac, dac_contents(j).name);
                    if exist(dir_prof,'dir') == 7
                        if verbose
                            fprintf('Deleting %s\n', dir_prof)
                        end
                        rmdir(dir_prof, 's')
                    end
                    % full snapshots contain meta, tech, and traj files
                    if ~keep
                        dir_float = sprintf('%s/%s/%s', snap_path_dac, ...
                            this_dac, dac_contents(j).name);
                        float_files = dir(dir_float);
                        fn_files = {float_files.name};
                        for f = 1:length(fn_files)
                            if endsWith(fn_files{f}, '.nc') && ...
                                    ~contains(fn_files{f}, 'prof')
                                delete([dir_float, '/', fn_files{f}])
                            end
                        end
                    end
                end
            end
        end
    end
end

% sort files into appropriate subdirectories - use the same setup
% as for files downloaded directly from the GDAC
if keep && ~strcmp(snap_type, 'bgc')
    subdirs = {'Index';'Profiles';'Meta';'Tech';'Traj'};
else
    subdirs = {'Index';'Profiles'};
end
for i = 1:length(subdirs)
    [status, message] = mkdir([snap_path, '/', subdirs{i}]);
    if ~status
        fprintf('%s/%s could not be created:\n%s', ...
            snap_path, subdirs{i}, message);
        return;
    end
end

% reorganize files into directories that are consistent with this toolbox
contents = dir(snap_path_dac);
for i = 1:length(contents)
    if contents(i).isdir
        if ~strcmp(contents(i).name, '.') && ...
                ~strcmp(contents(i).name, '..')
            if strcmp(snap_type, 'bgc')
                dac_files = sprintf('%s/%s/*prof.nc', snap_path_dac, ...
                    contents(i).name);
            else
                % each dac directory has subdirectories for each float
                % meta, tech, and traj files are included if they were kept
                dac_files = sprintf('%s/%s/*/*.nc', snap_path_dac, ...
                    contents(i).name);
            end
            all_files = dir(dac_files);
            for f = 1:length(all_files)
                % first check if file will be kept or deleted
                if ~isempty(dac) && strcmp(snap_type, 'bgc')
                    keep_dac = 0;
                    for d = 1:length(dac)
                        if endsWith(all_files(f).folder, dac{d})
                            keep_dac = 1;
                            break;
                        end
                    end
                    if ~keep_dac
                        delete([all_files(f).folder, '/', ...
                            all_files(f).name])
                        continue
                    end
                end
                if ~isempty(floats)
                    keep_this = 0;
                    for fl = 1:length(floats)
                        if startsWith(all_files(f).name, ...
                                sprintf('%d_', floats(fl)))
                            keep_this = 1;
                            break;
                        end
                    end
                    if ~keep_this
                        delete([all_files(f).folder, '/', ...
                            all_files(f).name])
                        continue;
                    end
                end
                if contains(all_files(f).name, 'prof')
                    new_path = sprintf('%s/Profiles', snap_path);
                elseif contains(all_files(f).name, 'meta')
                    new_path = sprintf('%s/Meta', snap_path);
                elseif contains(all_files(f).name, 'tech')
                    new_path = sprintf('%s/Tech', snap_path);
                elseif contains(all_files(f).name, 'traj')
                    new_path = sprintf('%s/Traj', snap_path);
                else
                    new_path = snap_path; % this should not happen
                end
                [status,message] = movefile([all_files(f).folder, '/', ...
                    all_files(f).name], new_path);
                if ~status
                    warning('An error occurred during an attempted move of "%s" to %s\n:%s\n', ...
                        [all_files(f).folder, '/', all_files(f).name], ...
                        snap_path, message)
                end
            end
        end
    else % regular file
        this_file = [contents(i).folder, '/', contents(i).name];
        if endsWith(contents(i).name, '.gz')
            gunzip(this_file)
            this_file = this_file(1:end-3); % strip out '.gz'
            % there is no return value for the gunzip operation
            if exist(this_file, 'file')
                delete([this_file, '.gz'])
            end
        end
        if contains(this_file, 'index')
            [status,message] = movefile(this_file, [snap_path, '/Index']);
        else
            [status,message] = movefile(this_file, snap_path);
        end
        if ~status
            warning('An error occurred during an attempted move of file "%s" to %s:\n%s\n', ...
                this_file, snap_path, message)
        end
    end
end % i

% the dac subdirectory should have no files left in it at this point
rmdir(snap_path_dac, 's')
