function good_variables = check_variables(variables, varargin)
% check_variables  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   good_variables = check_variables(variables, varargin)
%
% DESCRIPTION:
%   This function checks if the specified variable(s) are available.
%   The available variables are returned.
%
% INPUT:
%   variables : cell array of variable(s) (or one variable as a string)
%
% OPTIONAL INPUT:
%   'warning', warning : if warning is not empty, this string will be
%               displayed with each name of a variable that is not
%               available
%
% OUTPUT:
%   good_variables : cell array of available variable(s)
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

global Settings;

% set defaults
warn = [];

% parse optional arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'warning')
        warn = varargin{i+1};
    else
        warning('unknown option: %s', varargin{i});
    end
end

% make sure Settings is initialized
if isempty(Settings)
    initialize_argo();
end

if ischar(variables)
    variables = cellstr(variables);
end

good_variables = variables;

for i = 1:length(variables)
    if ~any(strcmp(variables{i}, Settings.avail_vars))
        if ~isempty(warn)
            warning('%s: %s', warn, variables{i});
        end
        good_variables(strcmp(variables{i}, good_variables)) = [];
    end
end
