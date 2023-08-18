function good_float_ids = show_timeseries(float_ids, variables, ...
    depth, varargin)
% show_timeseries  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   good_float_ids = show_timeseries(float_ids, variables, ...
%       depth, varargin)
%
% DESCRIPTION:
%   This an intermediary function that downloads and loads data for the
%   given float(s) and calls plot_timeseries to create the plot(s).
%
% INPUTS:
%   float_ids      : WMO ID(s) of the float(s); or the Data struct as 
%                    returned by the load_float_data function
%   variables      : cell array of variable(s) (i.e., sensor(s)) to show
%   depth          : array of depth levels to plot
%
% OPTIONAL INPUTS:
%   'end',end_date : end date (in one of the following formats:
%                   [YYYY MM DD HH MM SS] or [YYYY MM DD])
%   'float_profs',fp : cell array with per-float indices of the profiles to
%                   be shown, as returned by select_profiles
%                   (ignored if Data struct was passed in as first argument)
%   'legend',legend: legend (string) can be 'yes' to show legend along with
%                   plot (default) or 'no'
%   'per_float',per_float : show time series separately for each float (1)
%                   or all in one plot (0); default: 1
%   'png',basename: if basename is not empty, png files will be created
%                   for all plots; the file names will be
%                   <basename>_<WMOID>_<variable>_<depth>dbar.png if
%                   per_float is 1 or <basename>_<variable>_<depth>dbar.png
%                   if per_float is 0
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
%                   default: 'no' (i.e., plot adjusted data if available
%                   for all selected floats);
%                   'no_strict': plot only adjusted data, skip floats
%                   that have only raw data available
%   'start',start_date : start date (in one of the following formats:
%                   [YYYY MM DD HH MM SS] or [YYYY MM DD])
%   'time_label',label : use either years ('y'), months ('m'), or days ('d');
%                   default depends on length of time shown:
%                   'd' for up to 60 days, 'm' for up to 18 months,
%                   'y' otherwise
%   'title',title : title for the plot (default: "Depth: .. dbar"); an
%                   empty string ('') suppresses the title
%   'var2',variable: if variable is not empty, time series of this second
%                   variable will be plotted; if it is the same type as the
%                   first variable (e.g., DOXY2 compared to DOXY), it will
%                   be plotted using the same axes; otherwise, the right
%                   axis will be used for the second variable;
%                   this option can only be used with 'per_float',1
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

% assign empty arrays to all return values in case of early return
good_float_ids = [];

if nargin < 3 || ~isnumeric(depth)
    warning('Usage: show_timeseries(float_ids, variables, depth, varargin)')
    return
end

if isempty(float_ids)
    warning('no floats specified')
    return
end

% set defaults
float_profs = [];
basename = [];
var2 = [];
varargpass= {};

% parse optional arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'float_profs')
        float_profs = varargin{i+1};
    elseif strcmpi(varargin{i}, 'png')
        basename = varargin{i+1};
    elseif strcmp(varargin{i}, 'var2')
        var2 = check_variables(varargin{i+1}, 'warning', ...
            'unknown sensor will be ignored');
    else
        if strcmpi(varargin{i}, 'qc')
            if min(varargin{i+1}) < 0 || max(varargin{i+1}) > 9
                warning('only QC flags 0..9 are allowed!')
                continue; % don't add it to varargpass
            end
        end
        varargpass = [varargpass, varargin(i), varargin(i+1)];
    end
end

% convert requested variable to cell array if necessary and
% discard unknown variables
variables = check_variables(variables, 'warning', ...
    'unknown sensor will be ignored');

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
    Data = [];% if float profiles were specified, make sure that there are no empty
    % arrays; if so, disregard these floats
    if ~isempty(float_profs)
        no_profs = cellfun(@isempty, float_profs);
        if any(no_profs)
            warning('No profiles specified for float(s) %s', ...
                num2str(float_ids(no_profs)))
            float_ids(no_profs) = [];
        end
    end

    % download Sprof files if necessary
    good_float_ids = download_multi_floats(float_ids);
end

if isempty(good_float_ids)
    warning('no valid floats found')
else
    avail_vars = check_float_variables(good_float_ids, variables, ...
        'warning', 'Not available in all specified floats');
    var2 = check_float_variables(good_float_ids, var2, 'warning', ...
        'Not available in all specified floats');
    if isempty(Data)
        [Data, Mdata] = load_float_data(good_float_ids, [avail_vars; var2], ...
            float_profs);
    end
    plot_timeseries(Data, Mdata, avail_vars, depth, basename, ...
        'var2', var2, varargpass{:});
end
