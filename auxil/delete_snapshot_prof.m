function delete_snapshot_prof(snap_path_dac, fn_tar, keep, verbose)
% delete_snapshot_prof  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   delete_snapshot_prof(snap_path_dac, fn_tar, keep, verbose)
%
% DESCRIPTION:
%   This function deletes the individual profile files of a full
%   snapshot file. The calling function must check if keep is set to 3,
%   which indicates to keep these files.
%   Unless keep is non-zero, meta, tech, and traj files will be deleted
%   as well.
%
% INPUTS:
%   snap_path_dac: path to the 'dac' directory
%   fn_tar: name of the dac tarball (e.g., 'aoml_core.tar.gz')
%   keep: if 0, tech, meta, and traj files are deleted,
%         if > 0, they will be kept
%   verbose: if > 0, show which directories are being deleted
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

% extract the name of the dac of this tarball
match_dac = regexp(fn_tar, '[a-z]+_', 'match', 'once');
this_dac = match_dac(1:end-1); % without the trailing underscore
dac_contents = dir([snap_path_dac, '/', this_dac]);
for j = 1:length(dac_contents)
    if dac_contents(j).isdir && ...
            ~strcmp(dac_contents(j).name, '.') && ...
            ~strcmp(dac_contents(j).name, '..')
        dir_prof = sprintf('%s/%s/%s/profiles', snap_path_dac, ...
            this_dac, dac_contents(j).name);
        if exist(dir_prof,'dir') == 7
            if verbose > 1
                fprintf('Deleting %s\n', dir_prof)
            end
            rmdir(dir_prof, 's')
        end
        % full snapshots contain meta, tech, and traj files,
        % which will be kept only if requested
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
