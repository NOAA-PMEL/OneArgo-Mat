function [path_index, path_index_bgc] = unzip_index_files(verbose)
% unzip_index_files  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   [path_index, path_index_bgc] = unzip_index_files(verbose)
%
% DESCRIPTION:
%   This function unzips the Sprof index file and for full snapshots 
%   also the full index file for one snapshot of Argo data.
%   It deletes the gzipped files afterwards.
%   Settings.snap_dir and Settings.snap_path
%   values are used to determine the path.
%
% INPUT: 
%   verbose: if > 0, show information about what is being done
%
% OUTPUTS:
%   path_index: path to the Sprof index file in case of BGC snapshots
%               or to the prof index file in case of full snapshots
%   path_index_bgc: path to the Sprof index file
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

% all snapshot tarballs contain the Sprof index file
path_index = sprintf('%s%sdac/argo_synthetic-profile_index.txt.gz', ...
    Settings.snap_dir, Settings.snap_path);
if verbose
    fprintf('Unzipping %s\n', path_index)
end
gunzip(path_index)
delete(path_index)
path_index = path_index(1:end-3); % strip out trailing '.gz'
path_index_bgc = path_index;

% full snapshot tarballs also contain the prof index file
if ~contains(Settings.snap_path, 'Bgc')
    index_file = sprintf('%sdac/ar_index_global_prof.txt.gz', ...
        Settings.snap_path);
    path_index = sprintf('%s%s', Settings.snap_dir, index_file);
    if verbose
        fprintf('Unzipping %s\n', path_index)
    end
    gunzip(path_index)
    delete(path_index)
    path_index = path_index(1:end-3); % strip out trailing '.gz'
end
