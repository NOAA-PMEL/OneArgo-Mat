function plot_maps(Data, Mdata, variables, depths, sz, basename, varargin)
% plot_maps  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   plot_maps(Data, Mdata, variables, basename, varargin)
%
% DESCRIPTION:
%   This function plots maps of values for the specified variable(s)
%   and float(s) at the given depth(s).
%
% PREREQUISITE:
%   (S)prof file(s) for the specified float(s) must exist locally.
%
% INPUTS:
%   Data      : struct that must contain the PRES field and the given
%               variables (_ADJUSTED fields are used if available)
%   Mdata     : struct that must contain the WMO_ID field
%   variables : cell array with names of the measured fields (e.g., {'DOXY'})
%   depths    : array of the depths to be plotted
%   sz        : sz defines the size of plotted points
%   basename  : if not empty, create png files of all plots;
%               the file names will be <variable>_<depth>dbar.png             
%
% OPTIONAL INPUTS:
%   'caxis',[cmin cmax] : specify the minimum and maximum value to be
%                     shown (for all variables and depths)
%                     default: automatic scaling for each variable and depth
%                     separately
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
%                     available)
%                     'no_strict': plot only adjusted data, skip floats
%                     that have only raw data available
%
% OUTPUT: None.
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
    warning(['Usage: plot_maps(Data, Mdata, variables, depth, basename, ', ...
        '[, varargin])'])
    return;
end

% set defaults
cax = [];
depth = [];
raw = 'no'; % plot adjusted data by default
qc_flags = 0:9; % use all data
title_add = '';

% parse optional arguments
for f = 1:2:length(varargin)-1
    if strcmpi(varargin{f}, 'caxis')
        cax = varargin{f+1};
    elseif strcmpi(varargin{f}, 'raw')
        raw = varargin{f+1};
    elseif strcmpi(varargin{f}, 'qc')
        qc_flags = varargin{f+1};
    else
        warning('unknown option: %s', varargin{f});
    end
end

floats = fieldnames(Data);
nfloats = length(floats);
if ~nfloats
    warning('no floats found in Data structure')
    return
end

float_ids = fieldnames(Mdata);
nvars = length(variables);
ndepths = length(depths);
nplots = ndepths * nvars;

if nplots > Settings.max_plots
    warning(['too many plots requested - use fewer depths and/or variables', ...
        newline, 'or increase Settings.max_plots if possible'], ...
        warn_insert)
    return
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
                    'available;\nits values will not be shown'], ...
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
        end
        if has_adj == nfloats
            variables{v} = [variables{v}, '_ADJUSTED'];
        end
    end
end

[lon_lim, lat_lim, Data] = get_lon_lat_lims(Data);
% Set lat and lon limits
latlim = [lat_lim(1)-Settings.pad_lat, lat_lim(2)+Settings.pad_lat];
lonlim = [lon_lim(1)-Settings.pad_lon, lon_lim(2)+Settings.pad_lon];
% Adjust limits outside range to minimum and maximum limits
latlim(latlim < -90) = -90;
latlim(latlim >  90) =  90;
% use 0..360 range if all points are near the dateline
if isfield(Data.(floats{1}), 'ALT_LON')
    lonlim(lonlim < 0) = 0;
    lonlim(lonlim >  360) =  360;
    use_alt_lon = 1;
else % using a range of -180..180
    lonlim(lonlim < -180) = -180;
    lonlim(lonlim >  180) =  180;
    use_alt_lon = 0;
end

% vertical interpolation to depths with regular intervals
for f = 1:nfloats
    Datai.(floats{f}) = depth_interp(Data.(floats{f}), qc_flags, 'raw', raw);
end

for v = 1:nvars
    if ~any(show_var(v,:))
        fprintf('None of the specified floats have adjusted %s values, skipping the plot\n', ...
            variables{v})
        continue
    end
    [long_name, units] = get_var_name_units(variables{v});
    for d = 1:ndepths
        for f = 1:nfloats
            [mini,idz] = min(abs(Datai.(floats{f}).PRES(:,1) - depths(d)));
            if mini < 2 % interpolated depths are at 2 dbar intervals
                break
            end
        end
        if mini > 2
            warning('The requested depth level (%d dbar) was not found and will be skipped.', ...
                depths(d))
            continue;
        end
        f1 = figure;
        if isfield(Settings, 'colormap')
            cmap = colormap(Settings.colormap);
        else
            cmap = colormap;
        end
        geoaxes;
        hold on;
        % Set geographic limits for figure
        geolimits(latlim,lonlim)
        geobasemap grayland
        if ~isempty(cax)
            caxis(cax);
        end
        for f = 1:nfloats
            if ~show_var(v,f)
                continue;
            end
            if size(Datai.(floats{f}).(variables{v}), 1) < idz
                warning('Float %s did not reach %d dbar', floats{f}, ...
                    depths(d))
                continue
            end
            if use_alt_lon
                lon{f} = Data.(floats{f}).ALT_LON(1,:);
            else
                lon{f} = Data.(floats{f}).LONGITUDE(1,:);
            end
            good_pts = ~isnan(Datai.(floats{f}).(variables{v})(idz,:));
            if ~isempty(good_pts)
                if sum(good_pts) == 3
                    % when there are exactly 3 good points, Matlab
                    % interprets the values as an RGB triplet instead of
                    % variable values; run a loop over them instead
                    idx_good_pts = find(good_pts);
                    for i = 1:3
                        geoscatter(Data.(floats{f}).LATITUDE(1,idx_good_pts(i)), ...
                            lon{f}(idx_good_pts(i)), ...
                            sz,Datai.(floats{f}).(variables{v})(idz,idx_good_pts(i)), '.');
                    end

                else
                    geoscatter(Data.(floats{f}).LATITUDE(1,good_pts), lon{f}(good_pts), ...
                        sz,Datai.(floats{f}).(variables{v})(idz,good_pts), '.');
                end
            end
        end
        title(sprintf('%s %s at %d dbar%s', long_name, units, ...
            depths(d), title_added{v}), 'FontSize', 14);
        colorbar
        if ~isempty(basename)
            fn_png = sprintf('%s_%s_%ddbar.png', basename, variables{v}, depths(d));
            print(f1, '-dpng', fn_png);
        end
    end
end
