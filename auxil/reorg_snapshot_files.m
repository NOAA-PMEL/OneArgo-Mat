function reorg_snapshot_files(full_snap_path, orig_dac, floats, ...
    float_dacs, snap_type, keep)
% reorg_snapshot_files  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   reorg_snapshot_files(full_snap_path, orig_dac, floats, ...
%       float_dacs, snap_type, keep)
%
% DESCRIPTION:
%   This function reorganizes the files for one snapshot of Argo data
%   so that the directory tree is the same as for GDAC files.
%
% INPUTS: 
%   full_snap_path: relative or absolute path to the main directory of
%         the currently selected snapshot
%   orig_dac: cell array of DAC names whose float files will be kept
%         (empty to keep all float files unless floats are specified)
%   floats: list of WMOIDs of floats whose files will be kept (in addition
%         to those from orig_dac DACs, if specified); can be empty
%   float_dacs: list of DACs handling the specified floats; can be empty
%   snap_type: either 'all', 'phys', or 'bgc'
%   keep: if 0, tech, meta, and traj files will be deleted by this
%         function; they will be kept for keep > 0
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

snap_path_dac = [full_snap_path, 'dac'];

% reorganize files into directories that are consistent with this toolbox
contents = dir(snap_path_dac);
for i = 1:length(contents)
    if contents(i).isdir
        if ~strcmp(contents(i).name, '.') && ...
                ~strcmp(contents(i).name, '..')
            if strcmp(snap_type, 'bgc')
                if (~isempty(orig_dac) || ~isempty(floats)) && ...
                        ~any(strcmp(contents(i).name, orig_dac)) && ...
                        ~any(strcmp(contents(i).name, float_dacs))
                    rmdir([snap_path_dac, '/', contents(i).name], 's')
                    continue
                end
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
                if ~isempty(floats) && ...
                        ~any(strcmp(orig_dac, contents(i).name))
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
                    new_path = sprintf('%s/Profiles', full_snap_path);
                elseif contains(all_files(f).name, 'meta')
                    new_path = sprintf('%s/Meta', full_snap_path);
                elseif contains(all_files(f).name, 'tech')
                    new_path = sprintf('%s/Tech', full_snap_path);
                elseif contains(all_files(f).name, 'traj')
                    new_path = sprintf('%s/Traj', full_snap_path);
                else
                    new_path = full_snap_path; % this should not happen
                end
                if keep || contains(all_files(f).name, 'prof')
                    [status,message] = movefile([all_files(f).folder, '/', ...
                        all_files(f).name], new_path);
                    if ~status
                        warning('An error occurred during an attempted move of "%s" to %s\n:%s\n', ...
                            [all_files(f).folder, '/', all_files(f).name], ...
                            full_snap_path, message)
                    end
                end
            end
        end
    else % regular file
        this_file = [contents(i).folder, '/', contents(i).name];
        % if 'keep' is set to 2 or 3, original gzipped tarball files will be
        % kept; there is no need or reason for gunzipping them
        if endsWith(contents(i).name, '.gz') && ~contains(this_file, 'tar')
            gunzip(this_file)
            this_file = this_file(1:end-3); % strip out '.gz'
            % there is no return value for the gunzip operation
            if exist(this_file, 'file')
                delete([this_file, '.gz'])
            end
        end
        if contains(this_file, 'index')
            [status,message] = movefile(this_file, [full_snap_path, '/Index']);
        elseif ~contains(this_file, 'tar')
            [status,message] = movefile(this_file, full_snap_path);
        end
        if ~status
            warning('An error occurred during an attempted move of file "%s" to %s:\n%s\n', ...
                this_file, full_snap_path, message)
        end
    end
end
