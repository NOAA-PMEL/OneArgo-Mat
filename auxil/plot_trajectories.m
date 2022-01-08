function plot_trajectories(Data, color, title1, fn_png, float_ids, ...
    lines, lgnd, sz)
% plot_trajectories  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   plot_trajectories(Data, color, title1, fn_png, float_ids, ...
%       lines, lgnd, sz)
%
% DESCRIPTION:
%   This function plots the trajectories of one or more specified float(s).
%
% PREREQUISITE:
%   Sprof file(s) for the specified float(s) must exist locally.
%
% INPUTS:
%   Data  : struct that must contain LONGITUDE and LATITUDE fields
%   color : either 'multiple' (different colors for different floats),
%           or any standard Matlab color descriptor ('r', 'k', 'b',
%           'g' etc.; all trajectories will be plotted in the same color)
%           color can also be 'dac'; in this case, the trajectories
%           are colored by the DAC responsible for the floats
%   title1: title of the plot
%   fn_png: if not empty, create a png image of the plot with this file name
%   float_ids: WMO IDs of the floats to be plotted
%   lines : lines (string) can be 'yes' to connect float positions
%           with a line or 'no' (default)
%   lgnd  : lgnd (string) can be 'yes' to show legend along with
%           plot (default) or 'no'
%   sz    : sz defines the size of plotted points (default: 36)
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

global Settings Float;

if nargin < 8
    warning(['Usage: plot_trajectories(Data, color, title1, fn_png, ',...
        'float_ids, lines, lgnd, sz)'])
    return;
end

% Determine which floats have been imported
floats = fieldnames(Data);
nfloats = numel(floats);

[lon_lim, lat_lim, Data] = get_lon_lat_lims(Data);
% Set lat and lon limits
latlim = [lat_lim(1)-5 lat_lim(2)+5];
lonlim = [lon_lim(1)-5 lon_lim(2)+5];
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

f1 = figure;
if isfield(Settings, 'colormap')
    cmap = colormap(Settings.colormap);
else
    cmap = colormap;
end

% determine colormap for multiple floats
if strcmp(color, 'dac')
    dacs = unique(Float.dac(ismember(Float.wmoid, float_ids)));
    indx = round(linspace(1, size(cmap, 1), length(dacs)));
    cmap = colormap(cmap(indx,:));
    fidx = arrayfun(@(x) find(Float.wmoid==x, 1), float_ids);
elseif strcmp(color, 'multiple') && ~strcmp(Settings.mapping, 'native')
    if nfloats == 1
        color = 'r'; % use red for single floats
    else
        ncolor = size(cmap, 1);
    end
end
lon = cell(nfloats, 1);
for i = 1:nfloats
    if use_alt_lon
        lon{i} = Data.(floats{i}).ALT_LON(1,:);
    else
        lon{i} = Data.(floats{i}).LONGITUDE(1,:);
    end
end
% use "geoscatter" as default
if strcmp(Settings.mapping, 'native')
    geoaxes;
    hold on;
    % Set geographic limits for figure
    geolimits(latlim,lonlim)
    geobasemap grayland
    if strcmp(color, 'dac')
        % plot one point per DAC to create the legend
        for i = 1:length(dacs)
            idx1 = find(strcmp(Float.dac(fidx), dacs{i}), 1);
            if use_alt_lon
                lon1 = Data.(floats{idx1}).ALT_LON(1,1);
            else
                lon1 = Data.(floats{idx1}).LONGITUDE(1,1);
            end
            geoscatter(Data.(floats{idx1}).LATITUDE(1,1), lon1, ...
                sz, cmap(i, :), '.', 'DisplayName', dacs{i})
        end
        legend(dacs,'location','eastoutside','AutoUpdate','off')
    end
    % Plot float trajectories on map (points)
    for i = 1:nfloats
        if strcmp(color, 'multiple')
            geoscatter(Data.(floats{i}).LATITUDE(1,:), lon{i}, sz, '.');
        elseif strcmp(color, 'dac')
            this_dac = Float.dac(find(Float.wmoid == float_ids(i), 1));
            cidx = find(strcmp(dacs, this_dac), 1);
            geoscatter(Data.(floats{i}).LATITUDE(1,:), lon{i}, sz, ...
                repmat(cmap(cidx, :), size(lon,2), 1), '.')
        else
            geoscatter(Data.(floats{i}).LATITUDE(1,:), lon{i}, sz, ...
                [color, '.']);
        end
    end
    % add legend
    if strcmp(lgnd,'yes') && ~strcmp(color, 'dac')
        legend(floats,'location','eastoutside','AutoUpdate','off')
    end
    % reset color order
    set(gca,'ColorOrderIndex',1);
    % Plot lines on map
    if strcmp(lines,'yes')
        for i = 1:nfloats
            % mask out long jumps between opposite sides of the plot
            lon{i}(abs(diff(lon{i})) > 300) = nan;
            l = geoplot(Data.(floats{i}).LATITUDE(1,:), lon{i}, 'k', ...
                'linewidth', 1);
            % send line to back
            uistack(l,'bottom');
            clear l
        end
    end
elseif strcmp(Settings.mapping, 'm_map') % use "m_map" if indicated
    if diff(latlim) > 15
        proj = 'robinson';
    else
        proj = 'lambert';
    end
    m_proj(proj, 'lon', lonlim, 'lat', latlim);
    hold on;
    % Plot float trajectories on map (points)
    for i = 1:nfloats
        if strcmp(color, 'multiple') && nfloats > 1
            cidx = round(1 + (ncolor-1) * (i-1) / (nfloats - 1));
            m_scatter(lon{i}, Data.(floats{i}).LATITUDE(1,:), ...
                round(sz/6), cmap(cidx,:), 'filled');
        else
            m_scatter(lon{i}, Data.(floats{i}).LATITUDE(1,:), ...
                round(sz/6), color, 'filled');
        end
        m_grid('box','fancy');
        m_coast('patch', [0.7 0.7 0.7]);
    end
    % neither legend() nor m_legend() works for m_scatter plots
    % reset color order
    set(gca,'ColorOrderIndex',1);
    % Plot lines on map
    if strcmp(lines,'yes')
        for i = 1:nfloats
            % mask out long jumps between opposite sides of the plot
            lon{i}(abs(diff(lon{i})) > 300) = nan;
            l = m_plot(lon{i}, Data.(floats{i}).LATITUDE(1,:), ...
                'k', 'linewidth', 1);
            % send line to back
            uistack(l,'bottom');
            % send line up a number of positions equal to nfloats
            uistack(l,'up',nfloats);
            clear l
        end
    end
else % "plain" plot
    hold on;
    % Plot float trajectories on map (points)
    for i=1:nfloats
        if strcmp(color, 'multiple') && nfloats > 1
            cidx = round(1 + (ncolor-1) * (i-1) / (nfloats - 1));
            scatter(lon{i}, Data.(floats{i}).LATITUDE(1,:), ...
                round(sz/5), cmap(cidx,:), 'filled');
        else
            scatter(lon{i}, Data.(floats{i}).LATITUDE(1,:), ...
                round(sz/5), color, 'filled');
        end
    end
    box on
    xlabel('Longitude')
    ylabel('Latitude')
    % add legend
    if strcmp(lgnd,'yes') &&  ~strcmp(color, 'dac')
        legend(floats,'location','eastoutside','AutoUpdate','off')
    end
    % reset color order
    set(gca,'ColorOrderIndex',1);
    % Plot lines on map
    if strcmp(lines,'yes')
        for i = 1:nfloats
            % mask out long jumps between opposite sides of the plot
            lon{i}(abs(diff(lon{i})) > 300) = nan;
            l = plot(lon{i}, Data.(floats{i}).LATITUDE(1,:), ...
                'k', 'linewidth', 1);
            % send line to back
            uistack(l,'bottom');
            clear l
        end
    end
end
hold off
title(title1)
if ~isempty(fn_png)
    print(f1, '-dpng', fn_png);
end
