function plot_sections(Data, Mdata, variables, nvars, plot_isopyc, ...
    plot_mld, time_label, depth, raw, obs, basename, varargin)
% plot_sections  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   plot_sections(Data, Mdata, variables, nvars, plot_isopyc, ...
%                 plot_mld, time_label, depth, raw, obs, basename, varargin)
%
% DESCRIPTION:
%   This function plots sections of one or more specified float(s) for
%   the specified variable(s).
%
% PREREQUISITE:
%   Sprof file(s) for the specified float(s) must exist locally.
%
% INPUTS:
%   Data        : struct that must contain the PRES field and the given
%                 variables (_ADJUSTED fields are used if available)
%   Mdata       : struct that must contain the WMO_ID field
%   variables   : cell array with names of the measured fields (e.g., DOXY)
%   nvars       : only the first "nvars" variables from the Data field
%                 will be plotted (if plotting isopycnal lines and/or
%                 mixed layer depth is requested, TEMP and PSAL may have
%                 been added to the "variables" cell array)
%   plot_isopyc : if set to 1, isopycnal lines will be plotted at default
%                 values (24:27); specific sigma levels can be specified
%                 as well, e.g.:
%                 25
%                 [25.5, 26.3]
%                 25.5:0.1:26
%                 if set to 0, no isopycnal lines will be plotted
%   plot_mld    : if set to 1 or 2, mixed layer depth will be plotted,
%                 using either a temperature (1) or density (2) criterion
%   time_label  : either years ('y'), months ('m'), or days ('d')
%   depth       : minimum and maximum depths to plot (an empty array
%                 signals the plotting of all available depths)
%   raw         : if 'yes', always use raw values;
%                 if 'no', use adjusted variables if available for
%                 all selected floats;
%                 if 'no_strict', plot only adjusted data, skip floats
%                 that have only raw data available
%   obs         : if 'on', add dots at the depths of observations
%   basename    : if not empty, create png files of all plots;
%                 the file names will be <basename>_<variable>.png
%
% OPTIONAL INPUTS:
%   'end',end_date : end date (in one of the following formats:
%                 [YYYY MM DD HH MM SS] or [YYYY MM DD])
%   'qc',flags  : show only values with the given QC flags (array)
%                 0: no QC was performed;
%                 1: good data;
%                 2: probably good data;
%                 3: probably bad data that are
%                    potentially correctable;
%                 4: bad data;
%                 5: value changed;
%                 6,7: not used;
%                 8: estimated value;
%                 9: missing value
%                 default setting: 0:9 (all flags)
%                 See Table 7 in Bittig et al.:
%                 https://www.frontiersin.org/files/Articles/460352/fmars-06-00502-HTML-r1/image_m/fmars-06-00502-t007.jpg
%   'start',start_date : start date (in one of the following formats:
%                 [YYYY MM DD HH MM SS] or [YYYY MM DD])
%   'title_add',string : string will be added to the end of all titles;
%                 default: empty string
%
% OUTPUT: None
%
% AUTHORS:
%   J. Sharp and H. Frenzel (UW-CICOES), A. Fassbender (NOAA-PMEL), N. Buzby (UW)
%
% CITATION:
%   H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2025. OneArgo-Mat:
%   A MATLAB toolbox for accessing and visualizing Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.6588041
%
% LICENSE: oneargo_mat_license.m
%
% DATE: APRIL 16, 2025  (Version 1.1.0)

global Settings;

if nargin < 11
    warning(['Usage: plot_sections(Data, Mdata, variables, nvars, ', ...
        'plot_isopyc, plot_mld, time_label, depth, raw, obs, ', ...
        'basename [, varargin])']);
    return
end

% set defaults
qc_flags = []; % if not changed, actual defaults will be assigned below
start_date = [];
end_date = [];
title_add = [];

% parse optional arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'qc')
        qc_flags = varargin{i+1};
    elseif strcmp(varargin{i}, 'start')
        start_date = check_datenum(varargin{i+1});
    elseif strcmp(varargin{i}, 'end')
        end_date = check_datenum(varargin{i+1});
    elseif strcmp(varargin{i}, 'title_add')
        title_add = varargin{i+1};
    end
end

if isempty(qc_flags)
    qc_flags = 0:9; % use everything
end

floats = fieldnames(Data);
float_ids = fieldnames(Mdata);
nfloats = length(floats);
% note that Data may contain additional variables (e.g., TEMP and PSAL,
% which are needed to compute density, but plotting their sections was
% not requested - they will not be counted in nvars)
nplots = nfloats * nvars;
if nplots > Settings.max_plots
    warn_msg1 = sprintf('too many (%d) plots requested - use fewer profiles and/or variables\n', nplots);
    warn_msg2 = sprintf('or increase Settings.max_plots (%d) if possible', ...
        Settings.max_plots);
    warning([warn_msg1, warn_msg2])
    return
end

calc_dens = ~isequal(plot_isopyc, 0);

% default: assume adjusted values are available
show_var = ones(nvars, nfloats);

% unless 'raw' is specified, plot adjusted data
if strncmpi(raw,'y',1)
    [title_added{1:nvars}] = deal(' [raw values]');
else
    [title_added{1:nvars}] = deal(''); % default; may change later
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
                    ' showing raw value profiles instead'], ...
                    variables{v}, floats{f});
                title_added{v} = ' [raw values]';
            end
        end
        if has_adj == nfloats
            variables{v} = [variables{v} '_ADJUSTED'];
        end
    end
end

% vertical interpolation to depths with regular intervals
for f = 1:nfloats
    varargs = {};
    if plot_mld == 1 % temperature criterion
        varargs = {'calc_mld_temp', 1};
    elseif plot_mld == 2 % density criterion
        varargs = {'calc_mld_dens', 1};
    end
    Datai = depth_interp(Data.(floats{f}), qc_flags, ...
        'calc_dens', calc_dens, 'raw', raw, varargs{:});
    for v = 1:nvars
        if ~show_var(v,f)
            continue;
        end
        if ~isfield(Datai, variables{v})
            warning('Sensor %s not found for float %s', ...
                variables{v}, floats{f})
            continue;
        end
        if isempty(find(isfinite(Datai.PRES) & ...
                isfinite(Datai.(variables{v})), 1))
            warning('no valid data found for %s of float %s', ...
                variables{v}, floats{f})
            continue;
        end
        if size(Datai.PRES,2) == 1
            disp('at least 2 profiles are needed for section plots')
            warning('%s of float %s has only one profile', ...
                variables{v}, floats{f})
            continue;
        end
        [long_name, units] = get_var_name_units(variables{v});
        f1 = figure('Renderer', 'painters', 'Position', [10*f 10*f 800 400]);
        set(gca,'fontsize',16);
        pcolor(Datai.TIME, Datai.PRES, Datai.(variables{v}));
        shading('flat');
        hold on
        if strcmp(obs,'on')
            index = ~isnan(Data.(floats{f}).(variables{v}));
            scatter(Data.(floats{f}).TIME(index), Data.(floats{f}).PRES(index),...
                1,'.k','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);
        end
        if ~isequal(plot_isopyc, 0) && isfield(Datai, 'DENS')
            try
                dens = Datai.DENS_ADJUSTED;
            catch
                dens = Datai.DENS;
            end
            if isequal(plot_isopyc, 1)
                iso_levels = 24:27;
            elseif length(plot_isopyc) == 1
                iso_levels = [plot_isopyc plot_isopyc];
            else
                iso_levels = plot_isopyc;
            end
            [C,h]=contour(Datai.TIME,Datai.PRES,dens-1000,...
                iso_levels,'black','LineWidth',1);
            clabel(C,h,iso_levels);
        end
        if plot_mld == 1
            plot(Datai.TIME,Datai.MLD_TEMP(1,:),'k','LineWidth',2);
        elseif plot_mld == 2 && isfield(Datai, 'MLD_DENS')
            % some old core floats don't have PSAL and therefore
            % density cannot be computed
            plot(Datai.TIME,Datai.MLD_DENS(1,:),'k','LineWidth',2);
        end
        hold off
        mod_xaxis_time(Datai.TIME(1,1), Datai.TIME(1,end), ...
            start_date, end_date, time_label)
        if length(depth) == 2
            ylim(depth);
        end
        title(sprintf('Float %d: %s %s%s%s', ...
            Mdata.(float_ids{f}).WMO_NUMBER, long_name, units, ...
            title_added{v}, title_add), 'FontSize', 14);
        ylabel('Pressure (dbar)');
        set(gca,'Ydir','reverse');
        colorbar;
        caxis([min(min(Datai.(variables{v}))), ...
            max(max(Datai.(variables{v})))]);
        if ~isempty(basename)
            fn_png = sprintf('%s_%d_%s.png', basename, ...
                Mdata.(float_ids{f}).WMO_NUMBER, variables{v});
            print(f1, '-dpng', fn_png);
        end
    end
end
