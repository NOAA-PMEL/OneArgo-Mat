%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show_deep_floats.m
%
% This example script for the BGC-Argo MATLAB toolbox finds floats
% with oxygen sensors that reach depths of at least 3000 db.
% Oxygen and temperature profiles are plotted for the float that
% has the most deep profiles with oxygen data.
%
% AUTHORS:
%   H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL), N. Buzby (UW)
%
% CITATION:
%   H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2022. OneArgo-Mat:
%   A MATLAB toolbox for accessing and visualizing Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.6588042
%
% LICENSE: oneargo_mat_license.m
%
% DATE: JUNE 1, 2022  (Version 1.0.1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% search globally and over the whole available time span, find floats that
% have an oxygen sensor and reach at least 3000 m
depth = 3000;
[float_ids, float_profs] = select_profiles([],[],[],[],'depth',depth,...
    'sensor','DOXY');
fprintf('%d floats with BGC sensors reach at least 3000 db\n', ...
    length(float_ids))

show_trajectories(float_ids, 'position', 'last', 'color', 'multiple', ...
    'title', 'Most recent positions of deep BGC floats', 'legend', 'no');

sensors = list_sensors(float_ids);

%% load the data for further analysis
Data = load_float_data(float_ids, {'DOXY';'TEMP'}, float_profs);

%% filter out profiles that do not have valid oxygen values below 3000 m
float_names = fieldnames(Data);
for f = 1:length(float_names)
    good_o2 = isfinite(Data.(float_names{f}).DOXY);
    press = Data.(float_names{f}).PRES;
    press(~good_o2) = nan;
    is_deep = max(press) > depth;
    float_profs{f}(~is_deep) = [];
end
% remove floats without any matching profiles
no_fp = cellfun(@isempty, float_profs);
float_profs(no_fp) = [];
float_ids(no_fp) = [];

%% plot the DOXY and TEMP profiles for the float with the most profiles
if any(strcmp(sensors, 'DOXY')) && any(strcmp(sensors, 'TEMP'))
    num_float_profs = cellfun(@length, float_profs);
    [~, idx] = max(num_float_profs);
    show_profiles(float_ids(idx), 'DOXY', 'var2', 'TEMP')
end
