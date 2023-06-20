function good_float_ids = show_maps(float_ids, variables, depths, varargin)
% show_profiles  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   good_float_ids = show_maps(float_ids, variables, depths, varargin)
%
% DESCRIPTION:
%   This an intermediary function that downloads profile(s) for the given
%   float(s) and calls plot_maps to create the plot(s).
%
% INPUTS:
%   float_ids  : WMO ID(s) of the float(s); or the Data struct as 
%                returned by the load_float_data function
%   variables  : cell array of variable(s) or one variable as a string
%   depths     : depth levels to show (dbar)
%
% OPTIONAL INPUTS:
%   'caxis',[cmin cmax] : specify the minimum and maximum value to be
%                   shown (for all variables and depths);
%                   default: automatic scaling for each variable and depth
%                   separately
%   'float_profs',fp : cell array with per-float indices of the profiles to
%                   be shown, as returned by select_profiles
%                   (ignored if Data struct was passed in as first argument)
%   'png',basename: if basename is not empty, png files will be created
%                   for all plots; the file names will be
%                   <basename>_<variable>_<depth>dbar.png
%   'qc',flags    : show only values with the given QC flags (as an array)
%                   0: no QC was performed;
%                   1: good data;
%                   2: probably good data;
%                   3: probably bad data that are potentially correctable;
%                   4: bad data;
%                   5: value changed;
%                   6,7: not used;
%                   8: estimated value;
%                   9: missing value
%                   default setting: 0:9 (all flags)
%                   See Table 7 in Bittig et al.:
%                   https://www.frontiersin.org/files/Articles/460352/fmars-06-00502-HTML-r1/image_m/fmars-06-00502-t007.jpg
%   'raw',raw     : plot raw, i.e., unadjusted data if set to 'yes';
%                   default: 'no' (i.e., plot adjusted data if available)
%                   'no_strict': plot only adjusted data, skip floats
%                   that have only raw data available
%   'size',sz     : sz (positive integer) defines the size of plotted
%                   points (default: 100)
%
% OUTPUT:
%   good_float_ids : array of the float IDs whose Sprof files were
%                   successfully downloaded or existed already
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

% assign empty arrays to the return value in case of early return
good_float_ids = [];

if isempty(float_ids)
    warning('no floats specified')
    return
end

% set defaults
if nargin < 3
    warning('Usage: show_maps(float_ids, variables, depths, varargin)');
    return
end
float_profs = [];
basename = [];
sz = 100;
varargpass= {};

% parse optional arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'float_profs')
        float_profs = varargin{i+1};
    elseif strcmpi(varargin{i}, 'size')
        if round(varargin{i+1}) > 0
            sz = round(varargin{i+1});
        else
            warning('size must be a positive integer')
        end
    elseif strcmpi(varargin{i}, 'png')
        basename = varargin{i+1};
    else
        if strcmpi(varargin{i}, 'qc')
            if min(varargin{i+1}) < 0 || max(varargin{i+1}) > 9
                warning('only QC flags 0..9 are allowed!')
                continue; % don't add it to varargpass
            end
        end
        varargpass = [varargpass, varargin{i:i+1}];
    end
end

if ~strcmp(Settings.mapping, 'native')
    warning('This function is only implemented for the native Matlab geofunctions')
    return
end

% convert requested variable to cell array if necessary and
% discard unknown variables
variables = check_variables(variables, 'warning', ...
    'unknown sensor will be  ignored');

% check if alternate first argument (Data instead of float_ids) was used
if isstruct(float_ids)
    Data = float_ids;
    clear float_ids;
    % need to construct good_float_ids as a numerical array and
    % the Mdata struct with WMO_NUMBER entries
    Mdata = struct();
    str_floats = fieldnames(Data);
    nfloats = length(str_floats);
    good_float_ids = nan(nfloats, 1);
    for f = 1:nfloats
        good_float_ids(f) = str2double(str_floats{f}(2:end));
        Mdata.(str_floats{f}).WMO_NUMBER = good_float_ids(f);
    end
else
    Data = [];
    % if float profiles were specified, make sure that there are no empty
    % arrays; if so, disregard these floats
    if ~isempty(float_profs)
        no_profs = cellfun(@isempty, float_profs);
        if any(no_profs)
            warning('No profiles specified for float(s):');
            disp(float_ids(no_profs))
            float_ids(no_profs) = [];
            float_profs(no_profs) = [];
        end
    end

    % download prof and Sprof files if necessary
    good_float_ids = download_multi_floats(float_ids);
end

if isempty(good_float_ids)
    warning('no valid floats found')
else
    avail_vars = check_float_variables(good_float_ids, variables, ...
        'warning', 'Not available in all specified floats');
    if isempty(Data)
        [Data, Mdata] = load_float_data(good_float_ids, avail_vars, ...
            float_profs);
    end
    plot_maps(Data, Mdata, avail_vars, depths, sz, basename, varargpass{:});
end
