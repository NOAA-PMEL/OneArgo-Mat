function avail_variables = check_float_variables(float_ids, variables, ...
    varargin)
% check_float_variables  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   avail_variables = check_float_variables(float_ids, variables, varargin)
%
% DESCRIPTION:
%   This function checks if the specified variable(s) are available for
%   all of the specified floats.
%   The available variables are returned.
%
% INPUT:
%   float_ids : array with WMO ID(s) of the float(s)
%   variables : cell array of variable(s) (or one variable as a string)
%
% OPTIONAL INPUT:
%   'warning', warning : if warning is not empty, this string will be
%               displayed with each name of a variable that is not
%               available
%
% OUTPUT:
%   avail_variables : cell array of available variable(s)
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

global Settings Float;

% set defaults
avail_variables = [];
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

if isempty(float_ids)
    warning('No floats specified!');
    return
end

if ischar(variables)
    variables = cellstr(variables);
end

% download proif and Sprof files if necessary
good_float_ids = download_multi_floats(float_ids);

avail_variables = variables; % default assumption: all are available

for f = 1:length(good_float_ids)
    filename = sprintf('%s%s', Settings.prof_dir, ...
        Float.file_name{Float.wmoid == good_float_ids(f)});
    info = ncinfo(filename); % Read netcdf information
    these_vars = {info.Variables.('Name')};
    for i = 1:length(avail_variables)
        if ~any(strcmp(variables{i}, these_vars))
            if ~isempty(warn)
                warning('%s: %s', warn, variables{i});
            end
            avail_variables(strcmp(variables{i}, avail_variables)) = [];
        end
    end
end
