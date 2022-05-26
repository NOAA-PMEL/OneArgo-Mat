%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comp_doxy_doxy2.m
%
% This example script for the BGC-Argo MATLAB toolbox finds floats with
% at least two oxygen sensors and creates comparison plots between them.
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
% DATE: MAY 26, 2022  (Version 1.3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% search globally since January 1, 2018 for floats that have at least two
% oxygen sensors
[doxy2_floats,doxy2_fp] = select_profiles([],[],[2018,1,1],[],...
    'sensor',{'DOXY';'DOXY2'});

fprintf('%d floats with at least two oxygen sensors were found\n', ...
    length(doxy2_floats));

% show all sensors available on these floats
list_sensors(doxy2_floats);

% show all sensors available on the first of these floats
fprintf('\nSensors on the first float (%d):', doxy2_floats(1))
list_sensors(doxy2_floats(1));

%% show the locations of the floats
show_trajectories(doxy2_floats, 'color', 'multiple', 'title', ...
    'Floats with at least two oxygen sensors');

%% show time series for both oxygen sensors at three different depth levels,
% create png files
show_timeseries(doxy2_floats, 'DOXY', [10,100,200], 'var2', 'DOXY2', ...
    'float_profs', doxy2_fp, 'png', 'timeseries');

%% show comparison profiles for the first float - since both variables
% are of the same type (oxygen), only one x axis is shown
% The mean of the first sensor is shown in black with individual profiles
% (type 1) or standard deviation (type 2) shown in gray.
% The mean of the second sensor is shown in dark blue with individual profiles
% (type 1) or standard deviation (type 2) shown in light blue.
% These colors can be modified easily in initialize_argo by changing
% the values for Settings.color_var1_mean etc.

% type 1: all profiles with mean
show_profiles(doxy2_floats(1), 'DOXY','var2','DOXY2','png','doxy_doxy2_all');
% type 2: mean with standard deviation
show_profiles(doxy2_floats(1), 'DOXY','var2','DOXY2','method','mean','png','doxy_doxy2_mean');

%% show a correlation plot between the two sensors for the first float
Data = load_float_data(doxy2_floats(1), {'DOXY';'DOXY2'}, doxy2_fp(1));

float_names = fieldnames(Data);

f1 = figure();

scatter(Data.(float_names{1}).DOXY(:), Data.(float_names{1}).DOXY2(:), ...
    3, Data.(float_names{1}).PRES(:), 'filled');

hcb = colorbar('vert');
set(get(hcb, 'label'), 'String', 'Pressure (db)')

% make x and y axis look the same
xl = xlim;
yl = ylim;
min_val = min(xl(1), yl(1));
max_val = max(xl(2), yl(2));
xlim([min_val, max_val]);
ylim([min_val, max_val]);
if length(yticks) > length(xticks)
    xticks(yticks);
else
    yticks(xticks);
end

% line of identical values
hold on
l = plot(min_val:max_val, min_val:max_val, 'k-');
% send line to back
uistack(l,'bottom');
hold off
box on;

xlabel('Oxygen sensor 1 (\mumol kg^{-1})')
ylabel('Oxygen sensor 2 (\mumol kg^{-1})')
title(sprintf('Oygen sensor comparison for float %d', doxy2_floats(1)));
