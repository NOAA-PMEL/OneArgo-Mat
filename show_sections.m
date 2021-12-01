function good_float_ids = show_sections(float_ids, variables, varargin)
% show_sections  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
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
%                (if not set: Settings.demo_float is used)
%   variables  : cell array of variable(s) (i.e., sensor(s)) to show 
%                (if not set: {'DOXY'} (=O2) is used)
%
% OPTIONAL INPUTS:
%   'float_profs',fp   : per-float indices of the profiles to be shown,
%                        as returned by select_profiles
%   'isopyc',isopyc    : plot isopycnal lines if non-zero (default: 1=on)
%                        using a value of 1 will result in plotting
%                        isopycnal lines at default values (24:27);
%                        specific sigma levels can be specified as well, e.g.:
%                        'isopyc',25
%                        'isopyc',[25.5, 26.3]
%                        'isopyc',25.5:0.1:26
%                        if set to 0, no isopycnal lines will be plotted
%   'mld',mld          : plot mixed layer depth, using either a 
%                        temperature criterion (mld=1) or a density
%                        criterion (mld=2); default: 0=off
%   'time_label',label : use either years ('y', by default) or months ('m')
%   'max_depth',depth  : maximum depth to be plotted (default: all)
%   'raw',raw          : plot raw, i.e., unadjusted data if set to 'yes';
%                        default: 'no' (i.e., plot adjusted data if available)
%   'obs',obs          : if 'on', add dots at the depths of observations
%                        default: 'on'; use 'off' to turn off this behavior
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
%   'png',basename     : if basename is not empty, png files will be created
%                        for all plots, the file names will be
%                        <basename>_<variable>.png
%
% OUTPUT:
%   good_float_ids : array of the float IDs whose Sprof files were
%                    successfully downloaded or existed already
%
% AUTHORS: 
%   H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL), N. Buzby (UW),
%   J. Plant, T. Maurer, Y. Takeshita (MBARI), D. Nicholson (WHOI),
%   and A. Gray (UW)
%
% CITATION:
%   H. Frenzel*, J. Sharp*, A. Fassbender, N. Buzby, J. Plant, T. Maurer,
%   Y. Takeshita, D. Nicholson, A. Gray, 2021. BGC-Argo-Mat: A MATLAB
%   toolbox for accessing and visualizing Biogeochemical Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.4971318.
%   (*These authors contributed equally to the code.)
%
% LICENSE: bgc_argo_mat_license.m
%
% DATE: DECEMBER 1, 2021  (Version 1.1)

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
time_label = 'y';
max_depth = []; % used as flag: plot all available depths
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
    elseif strcmpi(varargin{i}, 'max_depth')
        max_depth = varargin{i+1};
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

% convert requested variable to cell array if necessary (string was used)
if ischar(variables)
    variables = cellstr(variables);
end

% download Sprof files if necessary
good_float_ids = download_multi_floats(float_ids);

if isempty(good_float_ids)
    warning('no valid floats found')
else
    nvars = length(variables);
    % add the necessary variables now, but don't plot their profiles
    if ~isequal(plot_isopyc, 0) || plot_mld
        if ~any(strcmp(variables,'TEMP'))
            variables{end+1} = 'TEMP';            
        end
        if ~any(strcmp(variables,'PSAL'))
            variables{end+1} = 'PSAL';            
        end
    end
    [Data, Mdata] = load_float_data(good_float_ids, variables, float_profs);    
    plot_sections(Data, Mdata, variables, nvars, plot_isopyc, plot_mld, ...
        time_label, max_depth, raw, obs, basename, varargpass{:})
end

