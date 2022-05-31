function sensors = list_sensors(float_ids, mode)
% list_sensors  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   sensors = list_sensors([float_ids] [, mode])
%
% DESCRIPTION:
%   This function lists all currently defined sensors by their
%   short names (e.g., DOXY) and long names (e.g., Oxygen) along with
%   the units. This can be done for all available floats or selected
%   floats.
%
% INPUT: None
%
% OPTIONAL INPUTS:
%   float_ids : array with WMO IDs of the floats to be considered
%               (default: all floats)
%   mode      : either 'all' (default; sensors that are available for all
%               profiles of the specified floats) or 'some' (sensors
%               that are available for some profiles of each of the
%               specified floats)
%
% OUTPUT:
%   sensors   : cell array with matching sensors
%
% AUTHORS:
%   H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL), N. Buzby (UW)
%
% CITATION:
%   H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2022. OneArgo-Mat:
%   A MATLAB toolbox for accessing and visualizing Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.6588042
%
% LICENSE: oneargo_mat_license.m
%
% DATE: JUNE 1, 2022  (Version 1.0.1)

global Float Settings;

% make sure Settings is initialized
if isempty(Settings)
    initialize_argo();
end

if nargin < 1
    % show all available sensors
    sensors = sort(Settings.avail_vars);
    fprintf('\nThese are the currently defined sensors:\n\n')
else
    if nargin < 2
        mode = 'all';
    end
    fidx = arrayfun(@(x) find(Float.wmoid==x, 1), float_ids, ...
        'UniformOutput', false);
    fidx(cellfun(@isempty, fidx)) = [];
    fidx=cell2mat(fidx);
    if isempty(fidx)
        warning('no valid float IDs specified')
        sensors = {};
        return
    end
    if strcmp(mode, 'all')
        sens = Float.min_sens;
    else
        sens = Float.max_sens;
    end
    vars = sens{fidx(1)};
    for f = 2:length(fidx)
        vars = intersect(vars, sens{fidx(f)});
    end
    sensors = sort(vars);
    fprintf(['\nThese are the sensors available on %s profiles of the ', ...
        'selected floats:\n\n'], mode);
end
nvars = length(sensors);
fprintf('%-20s   %-42s  %s\n', 'Short name', 'Long name', 'Units');
fprintf('%s\n', repelem('-', 80));
for i = 1:nvars
    [long_name, units] = get_var_name_units(sensors{i});
    fprintf('%-20s   %-42s  %s\n', sensors{i}, long_name, ...
        strrep(strrep(strrep(units, '{', ''), '}', ''), '\mu', 'u'))
end
