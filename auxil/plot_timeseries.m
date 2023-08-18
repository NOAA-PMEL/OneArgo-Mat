function plot_timeseries(Data, Mdata, variables, depth, basename, varargin)
% plot_timeseries  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   plot_timeseries(Data, Mdata, variables, depth, basename, varargin)
%
% DESCRIPTION:
%   This function plots time series of one or more specified float(s) for
%   the specified variable(s).
%
% PREREQUISITE:
%   Sprof file(s) for the specified float(s) must exist locally.
%
% INPUTS:
%   Data      : struct that must contain the PRES field and the given
%               variables (_ADJUSTED fields are used if available)
%   Mdata     : struct that must contain the WMO_ID field
%   variables : cell array with names of the measured fields (e.g., DOXY)
%   depth     : array of depth levels to plot
%   basename  : if not empty, create png files of all plots;
%               if per_float is used, the file names will be
%               <basename>_<WMOID>_<variable>.png,
%               if per_float is not used, the file names will be
%               <basename>_<variable>.png
%
% OPTIONAL INPUTS:
%   'end',end_date  : end date (in one of the following formats:
%                     [YYYY MM DD HH MM SS] or [YYYY MM DD])
%   'legend',legend : legend (string) can be 'yes' to show legend along with
%                     plot (default) or 'no'
%   'per_float',per_float : show time series separately for each float (1)
%                     or all in one plot (0); default: 1
%   'qc',flags      : show only values with the given QC flags (array)
%                     0: no QC was performed;
%                     1: good data;
%                     2: probably good data;
%                     3: probably bad data that are potentially correctable;
%                     4: bad data;
%                     5: value changed;
%                     6,7: not used;
%                     8: estimated value;
%                     9: missing value
%                     default setting: 0:9 (all flags)
%                     See Table 7 in Bittig et al.:
%                     https://www.frontiersin.org/files/Articles/460352/fmars-06-00502-HTML-r1/image_m/fmars-06-00502-t007.jpg
%   'raw',raw       : plot raw, i.e., unadjusted data if set to 'yes';
%                     default: 'no' (i.e., plot adjusted data if
%                     available for all selected floats);
%                     'no_strict': plot only adjusted data, skip floats
%                     that have only raw data available
%   'start',start_date : start date (in one of the following formats:
%                     [YYYY MM DD HH MM SS] or [YYYY MM DD])
%   'time_label',label : use either years ('y'), months ('m'), or days ('d');
%                     default depends on length of time shown:
%                     'd' for up to 60 days, 'm' for up to 18 months,
%                     'y' otherwise
%   'title',title   : title for the plot (default: "Depth: .. dbar"); an
%                     empty string ('') suppresses the title
%   'var2',variable : if variable is not empty, profiles of this second
%                     variable will be plotted; if it is the same type as the
%                     first variable (e.g., DOXY2 compared to DOXY), it will
%                     be plotted using the same axes; otherwise, the right
%                     axis will be used for the second variable
%
% OUTPUT: None
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

if nargin < 5
    warning(['Usage: plot_timeseries(Data, Mdata, variables, depth, ', ...
        'basename [, varargin])'])
    return;
end

% set defaults
per_float = 1; % show time series for each float in a separate plot
raw = 'no'; % plot adjusted data by default
title1 = []; % generate default titles
qc_flags = 0:9; % use all data
var2_orig = [];
lgnd = 'yes';
start_date = [];
end_date = [];
time_label = []; % used as flag

% parse optional arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'per_float')
        per_float = varargin{i+1};
    elseif strcmpi(varargin{i}, 'raw')
        raw = varargin{i+1};
    elseif strcmpi(varargin{i}, 'title')
        title1 = varargin{i+1};
    elseif strcmpi(varargin{i}, 'qc')
        qc_flags = varargin{i+1};
    elseif strcmpi(varargin{i}, 'var2')
        var2_orig = char(varargin{i+1});
    elseif strcmpi(varargin{i}, 'legend')
        lgnd = varargin{i+1};
    elseif strcmpi(varargin{i}, 'time_label')
        time_label = varargin{i+1};
    elseif strcmp(varargin{i}, 'start')
        start_date = check_datenum(varargin{i+1});
    elseif strcmp(varargin{i}, 'end')
        end_date = check_datenum(varargin{i+1});
    else
        warning('unknown option: %s', varargin{i});
    end
end

if ~isempty(var2_orig) && ~per_float
    warning('var2 can only be used with per float time series plots')
    return
end

floats = fieldnames(Data);
nfloats = length(floats);
if ~nfloats
    warning('no floats found in Data structure')
    return
end

float_ids = fieldnames(Mdata);
nvars = length(variables);
ndepths = length(depth);
if per_float
    nplots = nfloats * nvars * ndepths;
    warn_insert = 'floats, ';
else
    nplots = nvars * ndepths;
    warn_insert = '';
end
if nplots > Settings.max_plots
    warning(['too many plots requested (%d) - use fewer %sdepths and/or variables', ...
        newline, 'or increase Settings.max_plots (%d) if possible'], ...
        nplots, warn_insert, Settings.max_plots)
    return
end

% default: assume adjusted values are available
show_var = ones(nvars, nfloats); 
if isempty(var2_orig)
    show_var2 = zeros(1, nfloats);
else
    [var2{1:nvars}] = deal(var2_orig); % default; may change later
    show_var2 = ones(1, nfloats);
end

[title_added{1:nvars}] = deal(''); % default; may change later
% unless 'raw' is specified, plot adjusted data
if strncmpi(raw,'y',1)
    if ~ischar(title1) || ~isempty(title1) % not added if title1 == ''
        [title_added{1:nvars}] = deal(' [raw values]');
    end
else
    for v = 1:nvars
        % if all floats have adjusted values available for a variable,
        % they will be used instead of raw values
        has_adj = 0;
        for f = 1:nfloats
            % the "_ADJUSTED" variable usually exists, but it may not
            % be filled with actual values
            if isfield(Data.(floats{f}),[variables{v}, '_ADJUSTED']) && ...
                    sum(isfinite(Data.(floats{f}).([variables{v}, '_ADJUSTED'])(:)))
                has_adj = has_adj + 1;
            elseif strcmpi(raw, 'no_strict')
                warning(['adjusted values for %s for float %s are not available;\n',...
                    'this timeseries will not be shown'], variables{v}, ...
                    floats{f});
                show_var(v,f) = 0;
                has_adj = has_adj + 1; % needed for check below
            else
                warning(['adjusted values for %s for float %s are not available,',...
                    ' showing raw values for all floats instead'], ...
                    variables{v}, floats{f});
                if ~ischar(title1) || ~isempty(title1)
                    title_added{v} = [title_added{v}, ' [raw values]'];
                end
                break
            end
            if ~isempty(var2_orig)
                if isfield(Data.(floats{f}),[var2_orig, '_ADJUSTED']) && ...
                        sum(isfinite(Data.(floats{f}).([var2_orig, '_ADJUSTED'])(:)))
                    has_adj = has_adj + 1;
                elseif strcmpi(raw, 'no_strict')
                    warning(['adjusted values for %s for float %s are not available;\n',...
                        'this timeseries will not be shown'], var2_orig, ...
                        floats{f});
                    show_var2(f) = 0;
                    has_adj = has_adj + 1; % needed for check below
                else
                    warning(['adjusted values for %s for float %s ', ...
                        'are not available, showing raw values for ', ...
                        'all floats instead'], var2_orig, floats{f});
                    if ~ischar(title1) || ~isempty(title1)
                        title_added{v} = ' [raw values]';
                    end
                    break
                end
            end
        end
        if has_adj == nfloats * (1 + ~isempty(var2_orig))
            variables{v} = [variables{v}, '_ADJUSTED'];
            if ~isempty(var2_orig)
                var2{v} = [var2_orig, '_ADJUSTED'];
            end
        end
    end
end

if ~isempty(basename)
    [var2_insert{1:nvars}] = deal('');
    if ~isempty(var2_orig)
        for v = 1:nvars
            var2_insert{v} = sprintf('_%s', var2{v});
        end
    end
end

% vertical interpolation to depths with regular intervals
for f = 1:nfloats
    Datai.(floats{f}) = depth_interp(Data.(floats{f}), qc_flags, 'raw', raw);
end

% determine indices matching requested depths
ndepths = length(depth);
idx = cell(nfloats, ndepths);
press = cell(nfloats, ndepths); % actual pressure levels used
for d = 1:ndepths
    for f = 1:nfloats
        [mini, idx{f,d}] = min(abs(Datai.(floats{f}).PRES(:,1) - depth(d)));
        press{f,d} = Datai.(floats{f}).PRES(idx{f,d},1);
        if mini > Settings.depth_tol
            warning(['Closest depth level to requested %.1f dbar for ', ...
                'float %s is %.1f dbar'], depth(d), floats{f}, press{f,d});
        end
    end
end

% these variables will only be needed if there are floats with just
% one data point
f1 = figure(); % dummy figure, only used to get the color order
colors = get(gca,'ColorOrder');
close(f1);
n_colors = size(colors,1);

for v = 1:nvars
    if ~any(show_var(v,:)) && ~per_float
        fprintf('None of the specified floats have %s values, skipping the plot\n', ...
            variables{v})
        continue
    end
    if ~isempty(var2_orig)
        same_var_type = strncmp(variables{v}, var2{v}, ...
            length(variables{v}) - ...
            9 * endsWith(variables{v}, '_ADJUSTED'));
    end
    for d = 1:ndepths
        if ~per_float
            % one figure per variable for all floats; make wide plot
            f1 = figure('Position',[20 20 800 400]);
            % reset color order
            set(gca,'ColorOrderIndex',1);
            hold(gca, 'on')
            min_time = Inf;
            max_time = -Inf;
            for f = 1:nfloats
                if show_var(v,f)
                min_time = min(min_time, Datai.(floats{f}).TIME(1,1));
                max_time = max(max_time, Datai.(floats{f}).TIME(1,end));
            end
        end
        end
        for f = 1:nfloats
            if ~show_var(v,f)
                continue;
            end
            if per_float
                % one figure per variable for each float; make wide plot
                f1 = figure('Position',[20 20 800 400]);
                set(gca,'ColorOrderIndex',1);
                hold(gca, 'on')
            end
            % good values between missing values will not be connected by
            % a line, so plot them individually as points as well
            if length(Datai.(floats{f}).TIME(1,:)) == 1
                index  = get(gca,'ColorOrderIndex');
                if index > n_colors
                    index = 1;
                end
                color1 = colors(index,:);
                sc1 = scatter(Datai.(floats{f}).TIME(1,:), ...
                    Datai.(floats{f}).(variables{v})(idx{f,d},:), 7, ...
                    color1, 'filled');
            else
                plt1 = plot(Datai.(floats{f}).TIME(1,:), ...
                    Datai.(floats{f}).(variables{v})(idx{f,d},:), ...
                    'LineWidth', 2);
                color1 = get(plt1, 'color');
            end
            if find(isnan(Datai.(floats{f}).(variables{v})(idx{f,d},:)), 1)
                sc1 = scatter(Datai.(floats{f}).TIME(1,:), ...
                    Datai.(floats{f}).(variables{v})(idx{f,d},:), 7, ...
                    color1, 'filled');
                set(get(get(sc1, 'Annotation'), 'LegendInformation'), ...
                    'IconDisplayStyle', 'off');
            end
            [long_name, units] = get_var_name_units(variables{v});
            yl1 = ylabel([long_name, ' ', units]);
            if show_var2(f)
                if ~same_var_type
                    yyaxis right;
                end
                if length(Datai.(floats{f}).TIME(1,:)) == 1
                    index  = get(gca,'ColorOrderIndex');
                    if index > n_colors
                        index = 1;
                    end
                    color2 = colors(index,:);
                    sc2 = scatter(Datai.(floats{f}).TIME(1,:), ...
                        Datai.(floats{f}).(var2{v})(idx{f,d},:), 7, ...
                        color2, 'filled');
                else
                    plt2 = plot(Datai.(floats{f}).TIME(1,:), ...
                        Datai.(floats{f}).(var2{v})(idx{f,d},:), ...
                        'LineWidth', 2);
                    color2 = get(plt2, 'color');
                end
                if find(isnan(Datai.(floats{f}).(var2{v})(idx{f,d},:)), 1)
                    sc2 = scatter(Datai.(floats{f}).TIME(1,:), ...
                        Datai.(floats{f}).(var2{v})(idx{f,d},:), 7, ...
                        color2, 'filled');
                    set(get(get(sc2, 'Annotation'), 'LegendInformation'), ...
                        'IconDisplayStyle', 'off');
                end
                if ~same_var_type
                    [long_name, units] = get_var_name_units(var2_orig);
                    yl2 = ylabel([long_name, ' ', units]);
                    set(yl2, 'color', color2);
                    ax2 = gca;
                    ax2.YColor = color2;
                    yyaxis left;
                    ax1 = gca;
                    set(yl1, 'color', color1);
                    ax1.YColor = color1;
                end
            end
            if per_float
                hold off
                box on;
                mod_xaxis_time(Datai.(floats{f}).TIME(1,1), ...
                    Datai.(floats{f}).TIME(1,end), start_date, end_date, ...
                    time_label)
                if isempty(title1) && ~ischar(title1) % i.e., []
                    if strcmp(lgnd, 'no') || show_var2(f)
                        % add float number to title instead
                        title(sprintf('Depth: %d dbar (%s)%s', press{f,d}, ...
                            floats{f}, title_added{v}));
                    else % legend shows float number
                        title(sprintf('Depth: %d dbar%s', press{f,d}, ...
                            title_added{v}));
                    end
                else
                    title([title1, title_added{v}]);
                end
                if strcmp(lgnd,'yes')
                    if ~show_var2(f)
                        legend(float_ids{f},'location','eastoutside',...
                            'AutoUpdate','off');
                    else
                        legend({strrep(variables{v},'_','\_'); ...
                            strrep(var2{v},'_','\_')}, ...
                            'location', 'eastoutside', 'AutoUpdate', 'off');
                    end
                end
                if ~isempty(basename)
                    fn_png = sprintf('%s_%d_%s%s_%ddbar.png', basename, ...
                        Mdata.(float_ids{f}).WMO_NUMBER, variables{v}, ...
                        var2_insert{v}, press{f,d});
                    print(f1, '-dpng', fn_png);
                end
            end
        end
        if ~per_float
            hold off
            box on;
            mod_xaxis_time(min_time, max_time, start_date, ...
                end_date, time_label)
            if isempty(title1) && ~ischar(title1) % i.e., []
                if strcmp(lgnd, 'yes')
                    title(sprintf('Depth: %d dbar%s', press{f,d}, ...
                        title_added{v}));
                else
                    title(sprintf('Depth: %d dbar (%s)%s', press{f,d}, ...
                        floats{f}, title_added{v}));
                end
            else
                title([title1, title_added{v}]);
            end
            if strcmp(lgnd,'yes')
                legend(float_ids(show_var(v,:)==1),'location','eastoutside',...
                    'AutoUpdate','off');
            end
            if ~isempty(basename)
                fn_png = sprintf('%s_%s%s_%ddbar.png', basename, variables{v}, ...
                    var2_insert{v}, press{f,d});
                print(f1, '-dpng', fn_png);
            end
        end
    end
end
