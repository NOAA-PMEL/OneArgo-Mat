function download_snapshot(varargin)
% download_snapshot  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   download_snapshot(varargin)
%
% DESCRIPTION:
%   This function downloads and unpacks one snapshot of Argo data.
%   Settings.default_type and Settings.use_snapshot
%   values are used to determine which snapshot will be downloaded.
%
% INPUT: None
%
% OPTIONAL INPUTS  (key,value pairs):
%   'dac, dac: keep only files from the specified DAC(s), e.g., 'AOML'
%         or {'jma';'coriolis'}
%   'force', force: if 0 (default), do not download a snapshot file if
%         at least *some* of its contents exist locally already (NOTE:
%         the function does not check if specific files requested with
%         'dac' or 'floats' are present);
%         if 1, download snapshot in any case - this should be used if
%         additional floats or dac files need to be unpacked
%   'floats', floats: keep only files for the float(s) with the specified
%         WMO ID(s)
%         If both 'dac' and 'floats' are specified, all files needed for
%         both specifications are downloaded.
%         NOTE: If none of the specified floats are found for the specified
%         type of snapshot and 'dac' is not specified, a warning will be
%         issued and no files will be unpacked.
%   'keep', keep: if 0 (default), all files that are not used by this
%         toolbox are deleted; if >0, tech, meta, and traj files are kept
%         from full snapshots;
%         if >1, most files are kept ('aux' and 'geo' sub directories,
%         full tarball, intermediate tarballs in addition to files listed
%         for 'keep' set to 1)
%         if 3, individual profile netcdf files for full snapshots are
%         kept as well;
%         note that settings of 2 and 3 require a lot of disk space,
%         esp. for full snapshots
%         (If 'dac' and/or 'floats' options are used, non-matching float
%         profiles will be deleted even if keep is set to 2 or 3.)
%         If BGC snapshots are used, the only file affected by 'keep'
%         is the snapshot tarball file itself - kept if 'keep' is 2 or 3.
%   'verbose', verbose: if 1, show more information about what is being
%         done; if 0, only show errors and warnings; default: 1
%         if 2, show even more output about progress
%
% OUTPUT: None
%
% AUTHORS:
%   H. Frenzel and J. Sharp (UW-CICOES), A. Fassbender (NOAA-PMEL), N. Buzby (UW)
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

% make sure Settings is initialized
if isempty(Settings)
    initialize_argo();
end

% set default values
snap_type = Settings.default_type;
snap_date = Settings.use_snapshots;
dac = [];
floats = [];
force = 0;
keep = 0;
verbose = 1; % the untarring is slow, keep user informed

% parse optional arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'dac')
        dac = lower(varargin{i+1});
    elseif strcmpi(varargin{i}, 'floats')
        floats = lower(varargin{i+1});
    elseif strcmpi(varargin{i}, 'force')
        force = varargin{i+1};
    elseif strcmpi(varargin{i}, 'keep')
        keep = varargin{i+1};
    elseif strcmpi(varargin{i}, 'verbose')
        verbose = varargin{i+1};
    end
end

if ~snap_date
    fprintf('\nSettings.use_snapshots is set to 0, so no snapshot will be downloaded.\n')
    fprintf('If you want to switch to using snapshots, update the setting in\n')
    fprintf('initialize_argo.m, then run initialize_argo.\n\n')
    return
end

determine_snapshot();
if isempty(Settings.snap_path)
    return
end

full_snap_path = [Settings.snap_dir, Settings.snap_path];
if exist(full_snap_path, 'dir') && ~force
    if snap_date == 1
        fprintf('Most recent snapshot was downloaded already\n');
    else
        fprintf('Snapshot for %d was downloaded already\n', snap_date);
    end
    return
end

% always use cells for DAC(s); convert char array to cell
if ischar(dac)
    dac = cellstr(dac);
end

% check if the requested snapshot tarball exists already
download_snap = 1;
if exist(Settings.snap_file, 'file')
    % check if size of previously downloaded file matches the expected size
    file = dir(Settings.snap_file);
    file_size = file.bytes;
    if file_size == Settings.snap_size
        download_snap = 0;
    else
        fprintf('"%s" was not downloaded correctly before.\n', ...
            Settings.snap_file)
        fprintf('Expected file size: %d bytes\n', Settings.snap_size)
        fprintf('Actual file size:   %d bytes\n', file_size)        
    end
end
if download_snap
    fprintf('Starting download of "%s" now\n', Settings.snap_file)
    fprintf('Its file size is %.1f GB - this may take a while.\n', ...
        Settings.snap_size*1e-9)
    websave(Settings.snap_file, Settings.snap_url);
    if ~exist(Settings.snap_file, 'file')
        fprintf('"%s" could not be downloaded\n', Settings.snap_file)
        return
    end
    % check if size of downloaded file matches the expected size
    file = dir(Settings.snap_file);
    file_size = file.bytes;
    if file_size ~= Settings.snap_size
        fprintf('"%s" was not downloaded correctly.\n', Settings.snap_file)
        fprintf('Expected file size: %d bytes\n', Settings.snap_size)
        fprintf('Actual file size:   %d bytes\n', file_size)
        return
    end
end

try
    tic
    fprintf('Untarring the snapshot, this may take a few minutes... ')
    untar(Settings.snap_file, Settings.snap_dir);
    fprintf('done!\n')
    toc
catch
    fprintf('error! Aborting...\n')
    return
end

[path_index, path_index_bgc] = unzip_index_files(verbose);

[float_dacs, float_is_bgc] = determine_float_dacs(...
    floats, path_index, path_index_bgc);
if ~isempty(floats) && isempty(float_dacs)
    warning('None of the specified floats was found')
    if isempty(dac)
        warning('Snapshot will not be unpacked!')
        return
    end
end

% the original tarball is still present after untarring
% wait until here with the deletion in case wrong float IDs were used
% accidentally
if keep < 2
    delete(Settings.snap_file); % delete it to free up disk space
end

% If both 'dac' and 'floats' are specified, all files needed for
% both specifications are downloaded.
orig_dac = dac; % keep track if dac was specified as an option
dac = union(dac, float_dacs);

% determine which tarballs need to be unpacked
if ~strcmp(snap_type, 'bgc') && ~isempty(dac)
    if isempty(floats)
        need_dac_core = ones(length(dac), 1);
        need_dac_bgc = ones(length(dac), 1) .* (strcmp(snap_type, 'all') + 0);
    else
        need_dac_core = zeros(length(dac), 1);
        need_dac_bgc = zeros(length(dac), 1);
        for d = 1:length(dac)
            if any(strcmp(dac{d}, orig_dac))
                need_dac_core(d) = 1;
                need_dac_bgc(d) = 1; % not checked below if type is 'phys'
            else
                % find the specified floats from the current dac
                is_dac = find(strcmp(dac{d}, float_dacs));
                if any(float_is_bgc(is_dac))
                    need_dac_bgc(d) = 1;
                end
                if any(~float_is_bgc(is_dac)) || strcmp(snap_type, 'phys')
                    need_dac_core(d) = 1;
                end
            end
        end
    end
end

% Individual (S)prof files are sorted by dac subdirectory, which is
% a different approach from this toolbox's normal way of organizing files.
% Therefore, files will be reorganized.
snap_path_dac = [full_snap_path, 'dac'];

% In the case of full snapshots, there are core and bgc tarballs for
% each DAC, which must be unpacked first.
if ~strcmp(snap_type, 'bgc')
    % contents of aux and geo directory trees are not used by this toolbox
    if keep < 2
        rmdir([full_snap_path, '/aux'], 's');
        rmdir([full_snap_path, '/geo'], 's');
    end
    tarballs = dir([snap_path_dac, '/*.tar.gz']);
    fn_tar = {tarballs.name};
    for i = 1:length(fn_tar)
        if strcmp(snap_type, 'phys') && ...
                endsWith(fn_tar{i}, 'bgc.tar.gz')
            if keep < 2
                if verbose
                    fprintf('Deleting "%s"\n', fn_tar{i})
                end
                delete([snap_path_dac, '/', fn_tar{i}])
            end
            % if keep is set to 2 or 3, the bgc tarballs are kept, but not
            % unpacked
            continue
        end
        if ~isempty(dac) && keep < 2
            keep_tar = 0;
            for d = 1:length(dac)
                if startsWith(fn_tar{i}, dac{d})
                    if endsWith(fn_tar{i}, 'bgc.tar.gz')
                        keep_tar = need_dac_bgc(d);
                    elseif endsWith(fn_tar{i}, 'core.tar.gz')
                        keep_tar = need_dac_core(d);
                    end
                    break
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
            fprintf('Untarring "%s" (%.0f GB) may take a few minutes... ', ...
                fn_tar{i}, 1e-9*file_size);
        elseif verbose
            fprintf('Untarring "%s"... ', fn_tar{i})
        end
        try
            untar([snap_path_dac, '/', fn_tar{i}], snap_path_dac);
            fprintf('done!\n')
            if keep < 2 % delete tarball ASAP to free up disk space
                delete([snap_path_dac, '/', fn_tar{i}])
            end
        catch
            fprintf('Error during untar!\n')
        end
        % the full snapshot archives include individual profile files,
        % which are not used by this toolbox
        if keep < 3
            if verbose == 1
                fprintf('Deleting individual profile files from %s\n', ...
                    fn_tar{i})
            end
            delete_snapshot_prof(snap_path_dac, fn_tar{i}, keep, verbose)
        end
    end
end

% sort files into appropriate subdirectories - use the same setup
% as for files downloaded directly from the GDAC
% start by creating subdirectories
if keep && ~strcmp(snap_type, 'bgc')
    subdirs = {'Index';'Profiles';'Meta';'Tech';'Traj'};
else
    subdirs = {'Index';'Profiles'};
end
for i = 1:length(subdirs)
    [status, message] = mkdir([full_snap_path, '/', subdirs{i}]);
    if ~status
        fprintf('%s/%s could not be created:\n%s', ...
            full_snap_path, subdirs{i}, message);
        return;
    end
end

reorg_snapshot_files(full_snap_path, orig_dac, floats, float_dacs, ...
    snap_type, keep)

% if keep is 3, everything is kept
if keep <= 1 || strcmp(snap_type, 'bgc')
    % the dac subdirectory should have no files left in it at this point
    rmdir(snap_path_dac, 's')
elseif keep == 2
    % there are many empty directories below the 'dac' directory,
    % but also the tarball files that will be kept
    dir_tar = [full_snap_path, '/tarballs'];
    mkdir(dir_tar)
    for i = 1:length(fn_tar)
        movefile([snap_path_dac, '/', fn_tar{i}], dir_tar)
    end
    % the dac subdirectory should have no files left in it at this point
    rmdir(snap_path_dac, 's')
    % rename the directory with the tarball files to 'dac'
    movefile(dir_tar, snap_path_dac)
end