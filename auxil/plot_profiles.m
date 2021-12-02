function [mean_prof, std_prof, mean_pres] = plot_profiles(Data, Mdata, ...
    variables, basename, varargin)
% plot_profiles  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   [mean_prof, std_prof, mean_pres] = plot_profiles(Data, Mdata, ...
%       variables, varargin)
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
%   'method',method : either 'all' (all profiles from each float are
%                     shown in one plot per variable) or 'mean'
%                     (mean and standard deviation across profiles);
%                     default is 'all'
%   'per_float',per_float : show profiles separately for each float (1)
%                     or all in one plot (0); default: 1 --
%                     either option can be used with 'all' and
%                     'mean' methods
%   'obs',obs       : plot markers at depths of observations (1);
%                     default: 0 (=off)
%   'raw',raw       : plot raw, i.e., unadjusted data if set to 'yes';
%                     default: 'no' (i.e., plot adjusted data if 
%                     available)
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
%   'title_add',text: add the given text to the end of the title
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
%   J. Sharp, H. Frenzel, A. Fassbender (NOAA-PMEL), N. Buzby (UW),
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
method = 'all'; % show all profiles per variable in one plot
per_float = 1; % show profiles for each float in a separate plot
obs = 'off'; % don't show observation points on each profile by default
raw = 'no'; % plot adjusted data by default
title_add = ''; % nothing added to title
qc_flags = 0:9; % use all data

% parse optional arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'method')
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
nplots = nfloats * nvars;
if nplots > Settings.max_plots
    warning(['too many plots requested - use fewer profiles and/or variables', ...
        newline, 'or increase Settings.max_plots if possible'])
    return
end

% unless 'raw' is specified, plot adjusted data
if strncmpi(raw,'y',1)
    title_add = ' [raw values]';
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
            else
                warning(['adjusted values for %s for float %s are not available,',...
                    ' showing raw values for all profiles instead'], ...
                    variables{v}, floats{f});
                title_add = ' [raw values]';
            end
        end
        if has_adj == nfloats
            variables{v} = [variables{v} '_ADJUSTED'];
        end
    end
end

% vertical interpolation to depths with regular intervals
for f = 1:nfloats
    Datai.(floats{f}) = depth_interp(Data.(floats{f}), qc_flags);
end

if per_float
    for f = 1:nfloats
        mean_pres{f} = Datai.(floats{f}).PRES(:,1);
    end
end
for v = 1:nvars
    if ~per_float
        [mean_prof{v}, std_prof{v}, mean_pres] = ...
            get_multi_profile_mean(Datai, variables{v});
        this_mean_prof = mean_prof{v};
        this_std_prof = std_prof{v};
        this_mean_pres = mean_pres;
        f1 = figure; % one figure per variable for all floats
        hold on
    end
    for f = 1:nfloats
        try
            PRES = Data.(floats{f}).PRES_ADJUSTED;
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
            f1 = figure; % one figure per variable for each float
            hold on
        end
        if strcmp(method, 'all')
            for p = 1:nprofs
                idx = isfinite(Data.(floats{f}).(variables{v})(:,p)) & ...
                    isfinite(PRES(:,p)) & ...
                    ismember(Data.(floats{f}).([variables{v},'_QC'])(:,p),...
                    qc_flags);
                if sum(idx)
                    plot(Data.(floats{f}).(variables{v})(idx,p), ...
                        PRES(idx,p),'k-');
                    if strcmp(obs,'on')
                        scatter(Data.(floats{f}).(variables{v})(idx,p), ...
                            PRES(idx,p),'k.');
                    end
                else
                    warning(['no valid observations matching selected',...
                        ' QC flags found for profile %d of float %s'], ...
                        p, floats{f});
                end
            end
            plot(this_mean_prof,this_mean_pres,'r-','linewidth',2);
        else
            plot(this_mean_prof,this_mean_pres,'k-','linewidth',2);
            plot(this_mean_prof - this_std_prof,this_mean_pres,...
                'b-','linewidth',2);
            plot(this_mean_prof + this_std_prof,this_mean_pres,...
                'b-','linewidth',2);
        end
        set(gca,'Ydir','reverse');
        [long_name, units] = get_var_name_units(variables{v});
        xlabel([long_name, ' ', units])
        ylabel('Pressure (dbar)')
        if per_float
            hold off
            title(sprintf('Float %d %s', ...
                Mdata.(float_ids{f}).WMO_NUMBER, title_add));
            if ~isempty(basename)
                fn_png = sprintf('%s_%d_%s.png', basename, ...
                    Mdata.(float_ids{f}).WMO_NUMBER, variables{v});
                print(f1, '-dpng', fn_png);
            end
        end
    end
    if ~per_float
        hold off
        if nfloats < 4
            ttitle = 'Floats';
            for f = 1:nfloats
                ttitle = sprintf('%s %d', ttitle, ...
                    Mdata.(float_ids{f}).WMO_NUMBER);
            end
        else
            ttitle = 'All selected floats';
        end
        title([ttitle, title_add]);
        if ~isempty(basename)
            fn_png = sprintf('%s_%s.png', basename, variables{v});
            print(f1, '-dpng', fn_png);
        end
    end
end
