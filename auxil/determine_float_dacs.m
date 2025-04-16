function [float_dacs, float_is_bgc] = determine_float_dacs(...
    floats, path_index, path_index_bgc)
% determine_float_dacs  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   [float_dacs, float_is_bgc] = determine_float_dacs(...
%       floats, path_index, path_index_bgc)
%
% DESCRIPTION:
%   This function determines which DACs handle the specified floats and
%   whether or not they are BGC floats.
%   This function (instead of list_dacs) is needed for the use with 
%   snapshots, when the index files are not yet available in their
%   final locations.
%
% INPUTS:
%   floats: array of WMO IDs
%   path_index: path to the index file of all floats
%   path_index_bgc: path to the index file of BGC floats
%
% OUTPUTS:
%   float_dacs  : cell array of DACs
%   float_is_bgc: array of 0s and 1s
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

if isempty(floats)
    float_dacs = [];
    float_is_bgc = [];
    return
end

floats = unique(floats); % make sure that there are no duplicates
% need to use index file to determine the DACs handling these floats;
% need a custom algorithm because list_dacs uses the Float struct,
% which isn't created yet at this point
index_table = readtable(path_index, 'NumHeaderLines',8,'Delimiter',',');
split_path = split(index_table.file,'/');
match = regexp(split_path(:,2), '\d+', 'match');
float_ids = str2double([match{:}]');
index_table_bgc = readtable(path_index_bgc, 'NumHeaderLines',8,'Delimiter',',');
split_path_bgc = split(index_table_bgc.file,'/');
match_bgc = regexp(split_path_bgc(:,2), '\d+', 'match');
float_ids_bgc = str2double([match_bgc{:}]');
float_dacs = cell(length(floats), 1);
float_is_bgc = nan(length(floats), 1);
for f = 1:length(floats)
    all_float_dacs = split_path(ismember(float_ids, floats(f)), 1); % needed below
    if isempty(all_float_dacs)
        warning('float %d not found', floats(f))
    else
        float_dacs{f} = all_float_dacs{1};
    end
    float_is_bgc(f) = any(ismember(float_ids_bgc, floats(f)));
end
float_dacs(cellfun(@isempty, float_dacs)) = [];
