function good_float_ids = show_sections(float_ids, variables, varargin)
% show_sections  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   good_float_ids = show_sections(float_ids, variables, varargin)
%
% DESCRIPTION:
%   This an intermediary function that downloads profile(s) for the given
%   float(s) and calls plot_sections to create the plot(s).
%
% INPUTS:
%   float_ids  : WMO ID(s) of one or more floats
%                (if not set: Settings.demo_float is used); or 
%                the Data struct as returned by the load_float_data function
%   variables  : cell array of variable(s) (i.e., sensor(s)) to show
%                (if not set: {'DOXY'} (=O2) is used)
%
% OPTIONAL INPUTS:
%   'end',end_date     : end date (in one of the following formats:
%                        [YYYY MM DD HH MM SS] or [YYYY MM DD])
%   'depth',[min max]  : mininum and maximum depth to plot (default: all)
%   'float_profs',fp   : per-float indices of the profiles to be shown,
%                        as returned by select_profiles
%                        (ignored if Data struct was passed in as first argument)
%   'isopyc',isopyc    : plot isopycnal lines if non-zero (default: 1=on)
%                        using a value of 1 will result in plotting
%                        isopycnal lines at default values (24:27);
%                        specific sigma_0 levels can be specified as well, e.g.:
%                        'isopyc',25
%                        'isopyc',[25.5, 26.3]
%                        'isopyc',25.5:0.1:26
%                        if set to 0, no isopycnal lines will be plotted
%   'mld',mld          : plot mixed layer depth, using either a
%                        temperature criterion (mld=1) or a density
%                        criterion (mld=2); default: 0=off
%   'obs',obs          : if 'on', add dots at the depths of observations
%                        default: 'on'; use 'off' to turn off this behavior
%   'png',basename     : if basename is not empty, png files will be created
%                        for all plots, the file names will be
%                        <basename>_<variable>.png
%   'qc',flags         : show only values with the given QC flags (as an array)
%                        0: no QC was performed;
%                        1: good data;
%                        2: probably good data;
%                        3: probably bad data that are potentially correctable;
%                        4: bad data;
%                        5: value changed;
%                        6,7: not used;
%                        8: estimated value;
%                        9: missing value
%                        default setting: 0:9 (all flags)
%                        See Table 7 in Bittig et al.:
%                        https://www.frontiersin.org/files/Articles/460352/fmars-06-00502-HTML-r1/image_m/fmars-06-00502-t007.jpg
%   'raw',raw          : plot raw, i.e., unadjusted data if set to 'yes';
%                        default: 'no' (i.e., plot adjusted data if available
%                        for all selected floats);
%                        'no_strict': plot only adjusted data, skip floats
%                        that have only raw data available
%   'start',start_date : start date (in one of the following formats:
%                        [YYYY MM DD HH MM SS] or [YYYY MM DD])
%   'time_label',label : label can be years ('y'), months ('m'), or days ('d');
%                        default depends on length of time shown:
%                        'd' for up to 60 days, 'm' for up to 18 months,
%                        'y' otherwise
% OUTPUT:
%   good_float_ids : array of the float IDs whose Sprof files were
%                    successfully downloaded or existed already
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

if isempty(float_ids)
    warning('no floats specified')
    return
end

% set defaults
if nargin < 2
    variables = {'DOXY'};
end
if ~nargin
    float_ids = Settings.demo_float;
end
float_profs = [];
plot_isopyc = 1;
plot_mld = 0;
time_label = [];
depth = []; % used as flag: plot all available depths
raw = 'no'; % plot adjusted data by default
obs = 'on'; % plot observations on section by default
basename = [];
varargpass= {};

% parse optional arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'float_profs')
        float_profs = varargin{i+1};
    elseif strcmpi(varargin{i}, 'isopyc')
        plot_isopyc = varargin{i+1};
    elseif strcmpi(varargin{i}, 'mld')
        plot_mld = varargin{i+1};
    elseif strcmpi(varargin{i}, 'time_label')
        time_label = varargin{i+1};
    elseif strcmpi(varargin{i}, 'depth')
        depth = varargin{i+1};
    elseif strcmpi(varargin{i}, 'raw')
        raw = varargin{i+1};
    elseif (strcmpi(varargin{i}, 'obs'))
        obs = varargin{i+1};
    elseif strcmpi(varargin{i}, 'png')
        basename = varargin{i+1};
    else
        varargpass = [varargpass, varargin{i:i+1}];
        if strcmpi(varargin{i}, 'qc')
            if min(varargin{i+1}) < 0 || max(varargin{i+1}) > 9
                warning('only QC flags 0..9 are allowed!')
            end
        end
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
    Data = [];
    % download Sprof files if necessary
    good_float_ids = download_multi_floats(float_ids);
end

if isempty(good_float_ids)
    warning('no valid floats found')
else
    avail_vars = check_float_variables(good_float_ids, variables, ...
        'warning', 'Not available in all specified floats');
    nvars = length(avail_vars);
    % add the necessary variables now, but don't plot their profiles
    if ~isequal(plot_isopyc, 0) || plot_mld
        if ~any(strcmp(avail_vars,'TEMP'))
            avail_vars{end+1} = 'TEMP';
        end
        if ~any(strcmp(avail_vars,'PSAL'))
            avail_vars{end+1} = 'PSAL';
        end
    end
    if isempty(Data)
        [Data, Mdata] = load_float_data(good_float_ids, avail_vars, float_profs);
    end
    plot_sections(Data, Mdata, avail_vars, nvars, plot_isopyc, plot_mld, ...
        time_label, depth, raw, obs, basename, varargpass{:})
end
