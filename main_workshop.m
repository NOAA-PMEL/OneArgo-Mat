%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main_workshop.m
% Driver routine for the GO-BGC workshop Matlab tutorial
% June 28-30, 2021
% It uses the MATLAB toolbox for accessing BGC-Argo float data.
%
% Demonstrates the downloading of BGC-Argo float data with sample plots,
% a discussion of available data, quality control flags etc.
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

%% Close figures, clean up workspace, clear command window
close all; clear; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 0: Initialize
% This function defines standard settings and paths and creates Index
% and Profiles folders in your current path. It also downloads the Sprof 
% index file from the GDAC to your Index folder. The Sprof index is 
% referenced when downloading and subsetting float data based on user 
% specified criteria in other functions.
initialize_argo();
do_pause();

%% Examine global structures
% These global structures contain a variety of useful variables for
% downloading and manipulating float data. 'Sprof' contains fields
% with information for each profile, 'Float' contains fields with
% information for each float, 'Settings' contains settings to be used in
% the backgroud during plotting, etc. Variables in the global structures 
% can be altered within the initialize_argo.m file.
global Sprof Float Settings;

% Example: Look at the profile ID numbers and available sensors for the
% profiles that have been executed by new GO-BGC float #5906439.
float_idx = (Float.wmoid == 5906439); % index for float #5906439
prof_ids = Float.prof_idx1(float_idx):Float.prof_idx2(float_idx) % profile IDs for float #5906439
dates = datestr(datenum(Sprof.date(prof_ids), 'yyyymmddHHMMSS')) % dates of each profile from float #5906439
sensors = unique(Sprof.sens(prof_ids)) % sensors available for float #5906439
do_pause();

clear float_idx prof_ids dates sensors % clean up workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 1: SOCCOM float
% In this exercise, we download the NetCDF file for a Southern Ocean  
% BGC float, inspect its contents, show the trajectory, plot profiles
% for unadjusted and adjusted data, and show the effect of adjustments 
% made to the nitrate concentrations.

%% Download NetCDF file for float #5904183, a SOCCOM float with multiple seasons under ice
WMO = 5904859; 
success = download_float(WMO);

%% Display attributes, dimensions, and variables available in the NetCDF
ncdisp(['./Profiles/' num2str(WMO) '_Sprof.nc'])
do_pause();

%% Extract informational data from the NetCDF
S = ncinfo(['./Profiles/' num2str(WMO) '_Sprof.nc']);
{S.Variables.Name}' % show the available variables
do_pause();

%% We see that NITRATE is available, so load it (along with TEMP and PSAL) from the NetCDF
[data,mdata] = load_float_data(WMO,... % specify WMO number
    {'PSAL','TEMP','NITRATE'}); % specify variables
data.(['F' num2str(WMO)]) % show data that has been loaded into MATLAB
do_pause();

%% Show the trajectory of the downloaded float
show_trajectories(WMO);
do_pause();

%% Show all profiles for salinity and nitrate from the downloaded float
% this plots the raw, unadjusted data, and includes multiple profiles 
% compromised by biofouling that has affected the optics.
show_profiles(WMO, {'PSAL';'NITRATE'},'obs','on','raw','yes');
% this plots the adjusted data.
show_profiles(WMO, {'PSAL';'NITRATE'},'obs','on');
% this plots the adjusted, good (qc flag 1) and probably-good (qc flag 2) data.
show_profiles(WMO, {'PSAL';'NITRATE'},'obs','on','qc',[1 2]);
do_pause();

%% Show sections for nitrate
% this shows the raw, unadjusted data (pcolor plot)
% mixed layer depth is shown based on the temperature threshold
% (set the value to 2 after 'mld' to use the density threshold instead)
show_sections(WMO, {'NITRATE'},...
    'mld', 1,...  % tells the function to plot mixed layer depth using T
    'raw','yes'); % tells the function to plot raw data

show_sections(WMO, {'NITRATE'},...
    'mld', 2,...  % tells the function to plot mixed layer depth using rho
    'raw','no'); % tells the function to plot adjusted data (that is the
% default and could be left out in this call)

show_sections(WMO, {'NITRATE'}, 'mld', 2,...
    'qc',[1 2]); % tells the function to plot good and probably-good data
do_pause();

%% Clean up the workspace
clear data mdata S success WMO ans
clc; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 2: Ocean Station Papa floats
% In this exercise, we define a region in the Northeast Pacific along with
% a duration of time, and identify the float profiles matching that
% criteria. We show the trajectories of all the matching floats and plot
% profiles that match the criteria for one of the floats.

%% Set limits near Ocean Station Papa from 2008 to 2018
latlim=[45 60];
lonlim=[-150 -135];
t1=[2008 1 1];
t2=[2018 12 31];

%% select profiles based on those limits with specified sensor (NITRATE)
[OSP_floats,OSP_float_profs] = select_profiles(lonlim,latlim,t1,t2,...
    'sensor','NITRATE',... % this selects only floats with nitrate sensors
    'outside','both'); % All floats that cross into the time/space limits
                       % are identified from the Sprof index. The optional 
                       % 'outside' argument allows the user to specify
                       % whether to retain profiles from those floats that
                       % lie outside the space limits ('space'), time
                       % limits ('time'), both time and space limits 
                       % ('both'), or to exclude all profiles that fall 
                       % outside the limits ('none'). The default is 'none'.

% display the number of matching floats and profiles
disp(' ');
disp(['# of matching profiles: ' num2str(sum(cellfun('length',...
    OSP_float_profs)))]);
disp(['# of matching floats: ' num2str(length(OSP_floats))]);
disp(' ');

%% Show trajectories for the matching floats
% This function downloads the specified floats from the GDAC (unless the
% files have already been downloaded) and then loads the data for plotting.
% Adding the optional input pair 'color','multiple' will plot different
% floats in different colors
show_trajectories(OSP_floats,...
    'color','multiple'); % this plots different floats in different colors
do_pause();

%% Show profile plots for the first of these matching floats
% Case #1: all profiles from one float (1)
show_profiles(OSP_floats(1), {'PSAL';'DOXY'});
% Case #2: mean and standard deviation of all profiles from one float (1)
show_profiles(OSP_floats(1), {'PSAL';'DOXY'},...
    'method','mean'); % this tells the function to just plot the mean profile
do_pause();

%% clean up the workspace
clear data S success WMO ans OSP_floats OSP_profiles latlim lonlim t1 t2
clc;close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 3: Hawaii floats
% In this exercise, we define a region near Hawaii along with a duration of
% time. Again, we identify the float profiles matching those criteria, show
% their trajectories, plot all the matching profiles on one figure, and
% show sections for the unadjusted and adjusted values of salinity and
% dissolved oxygen.

%% Set limits near Hawaii from 2017 to 2019
latlim=[22 26];
lonlim=[-160 -155];
t1=[2017 1 1];
t2=[2019 12 31];

%% Select floats and profiles based on those limits
[HW_floats,HW_float_profs] = select_profiles(lonlim,latlim,t1,t2,...
    'outside','none'); % exclude profiles outside the time and space limits

% display the number of matching floats and profiles
disp(' ');
disp(['# of matching profiles: ' num2str(sum(cellfun('length',...
    HW_float_profs)))]);
disp(['# of matching floats: ' num2str(length(HW_floats))]);
disp(' ');

%% Show trajectories for the matching floats, along with the geo limits
% This function downloads the specified floats from the GDAC (unless the
% files have already been downloaded) and then loads the data for plotting.
% Adding the optional input pair 'color','multiple' will plot different
% floats in different colors
show_trajectories(HW_floats,'color','multiple');

% show domain of interest
hold on; 
if strcmp(Settings.mapping, 'native')
    geoplot([latlim(1) latlim(2) latlim(2) latlim(1) latlim(1)],...
        [lonlim(1) lonlim(1) lonlim(2) lonlim(2) lonlim(1)],...
        'k','linewidth',2);
elseif strcmp(Settings.mapping, 'm_map')
    m_plot([lonlim(1) lonlim(1) lonlim(2) lonlim(2) lonlim(1)],...
        [latlim(1) latlim(2) latlim(2) latlim(1) latlim(1)],...
        'k','linewidth',2);
else
    plot([lonlim(1) lonlim(1) lonlim(2) lonlim(2) lonlim(1)],...
        [latlim(1) latlim(2) latlim(2) latlim(1) latlim(1)],...
        'k','linewidth',2);
end
hold off;
do_pause();

%% Show trajectories for the matching profiles from each float, along with the geo limits
% Adding the optional input of 'float_profs' with the per-float profile numbers given by
% the select_profiles function will plot only the locations of those
% specified profiles from the specified floats
show_trajectories(HW_floats,'color','multiple','float_profs',HW_float_profs);

% show domain of interest
hold on; 
if strcmp(Settings.mapping, 'native')
    geoplot([latlim(1) latlim(2) latlim(2) latlim(1) latlim(1)],...
        [lonlim(1) lonlim(1) lonlim(2) lonlim(2) lonlim(1)],...
        'k','linewidth',2);
elseif strcmp(Settings.mapping, 'm_map')
    m_plot([lonlim(1) lonlim(1) lonlim(2) lonlim(2) lonlim(1)],...
        [latlim(1) latlim(2) latlim(2) latlim(1) latlim(1)],...
        'k','linewidth',2);
else
    plot([lonlim(1) lonlim(1) lonlim(2) lonlim(2) lonlim(1)],...
        [latlim(1) latlim(2) latlim(2) latlim(1) latlim(1)],...
        'k','linewidth',2);
end
hold off;
do_pause();

%% Show matching profiles from all floats
% show profiles (from all floats) within specified domain and times
show_profiles(HW_floats, {'PSAL';'DOXY'},'float_profs',HW_float_profs,...
    'per_float',0,...
    'qc',[1 2]);  % tells the function to plot good and probably-good data

do_pause();

%% Show only matching profiles from September
[~, ~, date] = get_lon_lat_time(HW_floats, HW_float_profs);

% determine profiles that occur in September for each float separately
for f = 1:length(HW_floats)
    dvec = datevec(date{f});
    month = dvec(:,2);
    HW_float_profs_Sep{f} = HW_float_profs{f}(month == 9);
end

show_profiles(HW_floats, {'PSAL';'DOXY'}, 'per_float', 0, ...
    'float_profs', HW_float_profs_Sep, ...
    'obs', 'on', ... % plot a marker at each observation
    'title_add', ' (September)',...  % add this to the title
    'qc', [1 2]); % apply QC flags
do_pause();

%% Show sections for pH and oxygen for the fifth float in the list of Hawaii floats
% this shows the raw, unadjusted data (pcolor plot)
% mixed layer depth is shown based on the temperature threshold
% (set the value to 2 after 'mld' to use the density threshold instead)
show_sections(HW_floats(5), {'PH_IN_SITU_TOTAL';'DOXY'},...
    'mld', 1,...   % tells the function to plot mixed layer depth using T
    'raw', 'yes'); % tells the function to plot raw (unadjusted) data
do_pause();

%% Show sections for pH and oxygen for the fifth float in the list of Hawaii floats
% this shows the adjusted data
show_sections(HW_floats(5), {'PH_IN_SITU_TOTAL';'DOXY'}, 'mld', 1,'raw', 'no');

%% clean up the workspace
clear all;
clc;close all
