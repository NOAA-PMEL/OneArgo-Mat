function plot_trajectories(Data, color)
% plot_trajectories  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   plot_trajectories(Data, color)
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
%
% AUTHORS: 
%   H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL),
%   J. Plant, T. Maurer, Y. Takeshita (MBARI), D. Nicholson (WHOI),
%   and A. Gray (UW)
%
% CITATION:
%   BGC-Argo-Mat: A MATLAB toolbox for accessing and visualizing
%   Biogeochemical Argo data,
%   H. Frenzel*, J. Sharp*, A. Fassbender, J. Plant, T. Maurer, 
%   Y. Takeshita, D. Nicholson, and A. Gray; 2021
%   (*These authors contributed equally to the code.)
%
% LICENSE: bgc_argo_mat_license.m
%
% DATE: June 15, 2021

global Settings;

if nargin < 2
    warning('Usage: plot_trajectories(Data, color)')
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
% use 0..360 range if all points are within 30 degrees of the dateline
if isfield(Data.(floats{1}), 'ALT_LON')
    lonlim(lonlim < 0) = 0;
    lonlim(lonlim >  360) =  360;
    use_alt_lon = 1;
else % using a range of -180..180
    lonlim(lonlim < -180) = -180;
    lonlim(lonlim >  180) =  180;
    use_alt_lon = 0;
end

figure;

% determine colormap for multiple floats
if strcmp(color, 'multiple') && ~strcmp(Settings.mapping, 'native')
    if nfloats == 1
        color = 'r'; % use red for single floats
    else
        cmap = colormap;
        ncolor = size(cmap, 1);
    end
end

if strcmp(Settings.mapping, 'native')
    geoaxes;
    hold on;
    % Set geographic limits for figure
    geolimits(latlim,lonlim)
    geobasemap grayland
    % Plot float trajectories on map
    for i=1:nfloats
        if use_alt_lon
            lon = Data.(floats{i}).ALT_LON(1,:);
        else
            lon = Data.(floats{i}).LONGITUDE(1,:);
        end
        if strcmp(color, 'multiple')
            geoscatter(Data.(floats{i}).LATITUDE(1,:), lon, '.');
        else
            geoscatter(Data.(floats{i}).LATITUDE(1,:), lon, ...
                [color, '.']);
        end
    end
    legend(floats,'location','eastoutside','AutoUpdate','off')
elseif strcmp(Settings.mapping, 'm_map')
    if diff(latlim) > 15
        proj = 'robinson';
    else
        proj = 'lambert';
    end
    m_proj(proj, 'lon', lonlim, 'lat', latlim);
    hold on;
    for i=1:nfloats
        if strcmp(color, 'multiple') && nfloats > 1
            cidx = round(1 + (ncolor-1) * (i-1) / (nfloats - 1));
            m_scatter(Data.(floats{i}).LONGITUDE(1,:),...
                Data.(floats{i}).LATITUDE(1,:),10,cmap(cidx,:),'filled');
        else
            m_scatter(Data.(floats{i}).LONGITUDE(1,:),...
                Data.(floats{i}).LATITUDE(1,:),10,color,'filled');
        end
        m_grid('box','fancy');
    end
    % neither legend() nor m_legend() works for m_scatter plots
else % "plain" plot
    hold on;
    for i=1:nfloats
        if strcmp(color, 'multiple') && nfloats > 1
            cidx = round(1 + (ncolor-1) * (i-1) / (nfloats - 1));
            scatter(Data.(floats{i}).LONGITUDE(1,:),...
                Data.(floats{i}).LATITUDE(1,:),10,cmap(cidx,:),'filled');
        else
            scatter(Data.(floats{i}).LONGITUDE(1,:),...
                Data.(floats{i}).LATITUDE(1,:),10,color,'filled');
        end
    end
    box on
    xlabel('Longitude')
    ylabel('Latitude')
    legend(floats,'location','eastoutside','AutoUpdate','off')
end
hold off
title('Float trajectories')
