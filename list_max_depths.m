function depths = list_max_depths(float_ids)
% list_max_depths  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   depths = list_max_depths(float_ids)
%
% DESCRIPTION:
%   This function lists the maximum depths reached by all profiles
%   and the deepest depth reached by one profile for each of the
%   selected floats. The values are also returned.
%
% INPUT:
%   float_ids : array with WMO IDs of the floats to be considered
%
% OUTPUT:
%   depths    : Nx3 array (N: number of floats); the first column contains
%               the maximum depths reached by all profiles per float,
%               the second column contains the deepest depth reached by
%               a profile per float;
%               the third column contains the float IDs.
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

if nargin < 1
    warning('Usage: list_max_depths(float_ids)')
    return
end

Data = load_float_data(float_ids, 'PRES');
floats = fieldnames(Data);
nfloats = length(floats);
if ~nfloats
    warning('no valid float numbers provided')
end

depths = nan(nfloats, 3);

fprintf('%-9s %-12s  %-12s%-6s\n', '  WMO ID', '  Max. depth', ...
    'Max. depth', 'Mode')
fprintf('%s(overall)   (all prof.) (pressure)\n', repelem(' ', 12));
fprintf('%s\n', repelem('-', 48));
for f = 1:nfloats
    this_float = strrep(char(floats{f}), 'F', '');
    depths(f, 3) = str2double(this_float);
    fprintf(' %-10s ', this_float);
    if isfield(Data.(floats{f}), 'PRES_ADJUSTED')
        mode = 'Adjusted';
        pres = Data.(floats{f}).PRES_ADJUSTED;
    else
        mode = 'Raw';
        pres = Data.(floats{f}).PRES;
    end
    max_prof = max(pres); % highest pressure for all profiles
    depths(f, 1) = max(max_prof);
    depths(f, 2) = min(max_prof);
    fprintf('%6.1f db  ', max(max_prof));
    fprintf(' %6.1f db ', min(max_prof));
    fprintf('  %s\n', mode)
end
