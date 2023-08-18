function [mean_prof, std_prof, mean_pres] = plot_profiles(Data, Mdata, ...
    variables, basename, varargin)
% plot_profiles  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   [mean_prof, std_prof, mean_pres] = plot_profiles(Data, Mdata, ...
%       variables, basename, varargin)
%
% DESCRIPTION:
%   This function plots profiles of one or more specified float(s) for
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
%   basename  : if not empty, create png files of all plots;
%               if per_float is used, the file names will be
%               <basename>_<WMOID>_<variable>.png,
%               if per_float is not used, the file names will be
%               <basename>_<variable>.png
%
% OPTIONAL INPUTS:
%   'depth',[min max] : minimum and maximum depth levels to plot
%   'method',method : either 'all' (all profiles from each float are
%                     shown in one plot per variable) or 'mean'
%                     (mean and standard deviation across profiles);
%                     default is 'all'
%   'obs',obs       : plot markers at depths of observations (1);
%                     default: 0 (=off)
%   'per_float',per_float : show profiles separately for each float (1)
%                     or all in one plot (0); default: 1 --
%                     either option can be used with 'all' and
%                     'mean' methods
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
%   'title_add',text: add the given text to the end of the title
%   'var2',variable : if variable is not empty, profiles of this second
%                     variable will be plotted; if it is the same type as the
%                     first variable (e.g., DOXY2 compared to DOXY), it will
%                     be plotted using the same axes; otherwise, right and
%                     top axes will be used for the second variable
%
% OUTPUTS:
%   mean_prof : mean across profiles (cell array of cell arrays of column
%               vectors ({variable}{float}) if per_float is set to 1,
%               cell array of column vectors ({variable}) if per_float
%               is set to 0)
%   std_prof  : standard deviation across profiles (same type as mean_prof)
%   mean_pres : mean pressure across profiles (cell array of column
%               vectors if per_float is set to 1,
%               column vector if per_float is 0)
%
% AUTHORS:
%   J. Sharp, H. Frenzel, A. Fassbender (NOAA-PMEL), N. Buzby (UW)
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

% empty return values in case of warnings
mean_prof = [];
std_prof = [];
mean_pres = [];

if nargin < 4
    warning(['Usage: plot_profiles(Data, Mdata, variables, basename, ', ...
        '[, varargin])'])
    return;
end

% set defaults
depth = [];
method = 'all'; % show all profiles per variable in one plot
per_float = 1; % show profiles for each float in a separate plot
obs = 'off'; % don't show observation points on each profile by default
raw = 'no'; % plot adjusted data by default
title_add = ''; % nothing added to title
qc_flags = 0:9; % use all data
var2_orig = [];

% parse optional arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'depth')
        if isnumeric(varargin{i+1}) && numel(varargin{i+1}) == 2
            depth = varargin{i+1};
            if depth(1) > depth(2) % ensure [min max] order
                depth = depth([2 1]);
            elseif depth(1) == depth(2)
                warning('min and max depth must be different values')
                depth = [];
            end
        else
            warning('min and max depth must be specified')
        end
    elseif strcmpi(varargin{i}, 'method')
        method = varargin{i+1};
    elseif strcmpi(varargin{i}, 'per_float')
        per_float = varargin{i+1};
    elseif strcmpi(varargin{i}, 'obs')
        obs = varargin{i+1};
    elseif strcmpi(varargin{i}, 'raw')
        raw = varargin{i+1};
    elseif strcmpi(varargin{i}, 'title_add')
        title_add = varargin{i+1};
    elseif strcmpi(varargin{i}, 'qc')
        qc_flags = varargin{i+1};
    elseif strcmpi(varargin{i}, 'var2')
        var2_orig = char(varargin{i+1});
    else
        warning('unknown option: %s', varargin{i});
    end
end

if ~strcmp(method, 'all') && ~strcmp(method, 'mean')
    warning('no such method! Only "all" and "mean" are possible.')
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
if per_float
    nplots = nfloats * nvars;
    warn_insert = 'floats and/or ';
else
    nplots = nvars;
    warn_insert = '';
end
if nplots > Settings.max_plots
    warning(['too many plots requested - use fewer %svariables', ...
        newline, 'or increase Settings.max_plots if possible'], ...
        warn_insert)
    return
end
if isempty(var2_orig)
    show_var2 = zeros(1, nfloats);
else
    [var2{1:nvars}] = deal(var2_orig); % default; may change later
    show_var2 = ones(1, nfloats);
end
% unless 'raw' is specified, plot adjusted data
show_var = ones(nvars, nfloats); % default: assume adjusted values are available
if strncmpi(raw,'y',1)
    [title_added{1:nvars}] = deal([title_add, ' [raw values]']);
else
    [title_added{1:nvars}] = deal(title_add);
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
                warning(['adjusted values for %s for float %s are not ',...
                    'available;\nits profiles will not be shown'], ...
                    variables{v}, floats{f});
                show_var(v,f) = 0;
                has_adj = has_adj + 1; % needed for check below
            else
                warning(['adjusted values for %s for float %s are not available,',...
                    ' showing raw values for all profiles instead'], ...
                    variables{v}, floats{f});
                title_added{v} = [title_added{v}, ' [raw values]'];
                break
            end
            if ~isempty(var2_orig)
                if isfield(Data.(floats{f}),[var2_orig, '_ADJUSTED']) && ...
                        sum(isfinite(Data.(floats{f}).([var2_orig, '_ADJUSTED'])(:)))
                    has_adj = has_adj + 1;
                elseif strcmpi(raw, 'no_strict')
                    warning(['adjusted values for %s for float %s are not ',...
                        'available;\nits profiles will not be shown'], ...
                        var2_orig, floats{f});
                    show_var2(f) = 0;
                    has_adj = has_adj + 1; % needed for check below
                else
                    warning(['adjusted values for %s for float %s ', ...
                        'are not available, showing raw values for ', ...
                        'all profiles instead'], var2_orig, floats{f});
                    title_added{v} = ' [raw values]';
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

if per_float
    mean_pres = cell(nfloats, 1);
    for f = 1:nfloats
        try
            mean_pres{f} = Datai.(floats{f}).PRES(:,1);
        catch % this happens when there are no valid values
            mean_pres{f} = [];
        end
    end
    mean_prof = cell(nvars, nfloats);
    std_prof = cell(nvars, nfloats);
else
    mean_prof = cell(nvars, 1);
    std_prof = cell(nvars, 1);
end
for v = 1:nvars
    if ~isempty(var2_orig)
        same_var_type = strncmp(variables{v}, var2{v}, ...
            length(variables{v}) - ...
            9 * endsWith(variables{v}, '_ADJUSTED'));
        if ~same_var_type && verLessThan('matlab', '9.7')
            warning(['Plotting profiles for two different variables ', ...
                'requires R2019b or later.']);
            continue;
        end
    end
    if ~per_float
        if ~any(show_var(v,:))
            fprintf('None of the specified floats have %s values, skipping the plot\n', ...
                variables{v})
            continue
        end
        [mean_prof{v}, std_prof{v}, mean_pres] = ...
            get_multi_profile_mean(Datai, variables{v});
        this_mean_prof = mean_prof{v};
        this_std_prof = std_prof{v};
        this_mean_pres = mean_pres;
        if ~isempty(var2_orig)
            [mean_prof_var2,std_prof_var2] = ...
                get_multi_profile_mean(Datai, var2{v});
        end
        f1 = figure; % one figure per variable for all floats
        if any(show_var2) && ~same_var_type
            [ax1, ax2] = create_tiled_layout();
        else
            ax1 = axes(f1);
            ax2 = ax1; % only used if ~isempty(var2_orig) && same_var_type
        end
        hold(ax1, 'on')
    end
    for f = 1:nfloats
        if ~show_var(v,f)
            continue;
        end
        try
            assert(strcmp(raw, 'no')) % use PRES if raw values are used
            PRES = Data.(floats{f}).PRES_ADJUSTED;
            good_vals = sum(isfinite(PRES));
            % a somewhat arbitrary criterion: at least half of the profiles
            % must have valid adjusted pressure values, or switch to
            % raw pressure
            assert(sum(good_vals) > 0.5 * length(good_vals));
        catch
            PRES = Data.(floats{f}).PRES;
        end
        nprofs = size(Data.(floats{f}).PRES,2);
        if per_float
            mean_prof{v}{f} = mean(Datai.(floats{f}).(variables{v}),2,...
                'omitnan');
            std_prof{v}{f} = std(Datai.(floats{f}).(variables{v}),[],2,...
                'omitnan');
            this_mean_prof = mean_prof{v}{f};
            this_std_prof = std_prof{v}{f};
            this_mean_pres = mean_pres{f};
            if ~isempty(var2_orig)
                mean_prof_var2 = mean(Datai.(floats{f}).(var2{v}),2,...
                    'omitnan');
                std_prof_var2 = std(Datai.(floats{f}).(var2{v}),[],2,...
                    'omitnan');
            end
            f1 = figure; % one figure per variable for each float
            if show_var2(f) && ~same_var_type
                [ax1, ax2] = create_tiled_layout();
            else
                ax1 = axes(f1);
                ax2 = ax1; % only used if ~isempty(var2_orig) && same_var_type
            end
            hold(ax1, 'on')
        end
        if strcmp(method, 'all')
            for p = 1:nprofs
                if ~plot_one_profile(ax1, ...
                        Data.(floats{f}).(variables{v})(:,p), ...
                        PRES(:,p), ...
                        Data.(floats{f}).([variables{v},'_QC'])(:,p), ...
                        qc_flags, strcmp(obs,'on'), Settings.color_var1_range)
                    warning(['no valid observations of %s matching selected',...
                        ' QC flags found for profile %d of float %s'], ...
                        variables{v}, p, floats{f});
                end
                if show_var2(f)
                    if ~plot_one_profile(ax2, Data.(floats{f}).(var2{v})(:,p), ...
                            PRES(:,p), ...
                            Data.(floats{f}).([var2{v},'_QC'])(:,p), ...
                            qc_flags, strcmp(obs,'on'), Settings.color_var2_range)
                        warning(['no valid observations of %s matching selected',...
                            ' QC flags found for profile %d of float %s'], ...
                            var2{v}, p, floats{f});
                    end
                end
            end
            plot(ax1,this_mean_prof,this_mean_pres,...
                'color',Settings.color_var1_mean,'linewidth',2);
            if show_var2(f)
                plot(ax2,mean_prof_var2,this_mean_pres,...
                    'color',Settings.color_var2_mean,'linewidth',2);
            end
        else
            plot(ax1,this_mean_prof,this_mean_pres,...
                'color',Settings.color_var1_mean,'linewidth',2);
            plot(ax1,this_mean_prof - this_std_prof,this_mean_pres,...
                'color',Settings.color_var1_range,'linewidth',2);
            plot(ax1,this_mean_prof + this_std_prof,this_mean_pres,...
                'color',Settings.color_var1_range,'linewidth',2);
            if show_var2(f)
                plot(ax2,mean_prof_var2,this_mean_pres,...
                    'color',Settings.color_var2_mean,'linewidth',2);
                plot(ax2,mean_prof_var2 - std_prof_var2,this_mean_pres,...
                    'color',Settings.color_var2_range,'linewidth',2);
                plot(ax2,mean_prof_var2 + std_prof_var2,this_mean_pres,...
                    'color',Settings.color_var2_range,'linewidth',2);
            end
        end
        set(ax1,'Ydir','reverse');
        [long_name, units] = get_var_name_units(variables{v});
        xlabel(ax1, [long_name, ' ', units])
        ylabel(ax1,'Pressure (dbar)')
        if ~isempty(depth)
            ylim(ax1, depth);
        end
        if show_var2(f) && ~same_var_type
            set(ax2, 'Ydir', 'reverse');
            [long_name, units] = get_var_name_units(var2_orig);
            xlabel(ax2, [long_name, ' ', units])
            ylabel(ax2, 'Pressure (dbar)')
            if isempty(depth)
                % make sure that both axes use the same pressure range
                yl1 = ylim(ax1);
                yl2 = ylim(ax2);
                ylim(ax1, [min(yl1(1), yl2(1)), max(yl1(2), yl2(2))])
                ylim(ax2, [min(yl1(1), yl2(1)), max(yl1(2), yl2(2))])
            else
                ylim(ax2, depth);
            end
        end
        if per_float
            hold off
            if ~show_var2(f) || same_var_type
                box(ax1, 'on');
            end
            title(sprintf('Float %d %s', ...
                Mdata.(float_ids{f}).WMO_NUMBER, title_added{v}));
            if ~isempty(basename)
                fn_png = sprintf('%s_%d_%s%s.png', basename, ...
                    Mdata.(float_ids{f}).WMO_NUMBER, variables{v}, ...
                    var2_insert{v});
                print(f1, '-dpng', fn_png);
            end
        end
    end
    if ~per_float
        hold off
        if ~show_var2(f) || same_var_type
            box(ax1, 'on')
        end
        if nfloats < 4
            ttitle = 'Floats';
            for f = 1:nfloats
                if show_var(v,f)
                    ttitle = sprintf('%s %d', ttitle, ...
                    Mdata.(float_ids{f}).WMO_NUMBER);
                end
            end
        else
            ttitle = 'All selected floats';
        end
        title([ttitle, title_added{v}]);
        if ~isempty(basename)
            fn_png = sprintf('%s%s.png', basename, variables{v}, ...
                var2_insert{v});
            print(f1, '-dpng', fn_png);
        end
    end
end
