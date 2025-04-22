%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FAIR_Tutorial.m
% Driver routine for the MATLAB toolbox for accessing One Argo float data.
%
% Initially written for the BGC Argo Access/Visualization Tutorial
% presented at FAIR Data Practices for Marine Ecological Time Series
% Workshop hosted April 22-25, 2025 at the Bermuda Institute of Ocean
% Sciences (BIOS).
%
% Inspired by the tutorial for the GO-BGC workshop presented by Alison
% Gray during June 28-30, 2021
%
% This tutorial demonstrates the accessibility of One Argo float data as
% well as the utilization of the BGC synthetic profiles (Sprofs) through
% exploration of meta data, sample profiles of BGC variables, and functions
% designed for BGC-Argo data visualization.
%
% This tutorial utilizes Sprof downloads from the GDAC.  It is recommended
% to use download_snapshot.m to obtain the monthly snapshot of Sprofs
% maintained at https://doi.org/10.17882/42182 for the purpose of research.
%
% Download_Snapshot_BGC_BATS.m provides an example of a snapshot download
% that pairs with this tutorial to demonstrate how to obtain Sprofs from
% the monthly snapshots. This method is not demonstrated in detail during
% this tutorial due to the length of time it takes to download the
% snapshots.
%
%
% TUTORIAL AUTHOR:
%   N. Guisewhite (MBARI)
%
% OneArgo-Mat AUTHORS:
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



%% Close figures, clean up workspace, clear command window
close all; clear; clear -global; clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 0.1: Initialize
% This function defines standard settings and paths and creates Index
% and Profiles folders in your current path. It also downloads the Sprof
% index file from the GDAC to your Index folder. The Sprof index is
% referenced when downloading and subsetting float data based on user
% specified criteria in other functions.
initialize_argo();
do_pause();


%% Exercise 0.2: Declare global structure
% The global structure 'Settings' contains settings to be used in
% the backgroud during plotting, etc. Its variables
% can be altered within the initialize_argo.m file.
global Settings;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 1.1: Look at the profile ID numbers and available sensors
% Let's look at profiles that have been executed by GO-BGC float #7901106.
WMO = 7901106;

% Load some data from this float (lon/lat, time cycle number etc.)
data = load_float_data(WMO);
str_floatnum = sprintf('F%d', WMO);

%% Exercise 1.2: Get dates of each profile from float #7901106
fprintf('\nA look at float %d - data are available for these dates:\n', WMO)

% Variables are stored in 2D arrays (number of levels x number of profiles)
% show the dates of the profiles (use the surface level; the time is the
% same for all pressure levels)
disp(datestr(data.(['F' num2str(WMO)]).TIME(1,:)))


%% Exercise 1.3: List sensors for 7901106
list_sensors(WMO); 
% This function simplifies the accessing of meta
% information stored in the netcdf, but you can still
% open and read the netcdf if you'd like.  That code is
% in the Extra Exercises section at the bottom of this
% script.


%% End of Exercise 1: Clean up the workspace
fprintf(['\nYou have finished Exercise 1! \nDuring the breakout portion, you can ',...
    'go back to this section and practice \nwith different WMOs.'...
    '\n\nProceeding will clear some variables and ',...
    'close any open plots.\n\n'])
do_pause();

clear WMO data str_floatnum ans % clean up workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 2.1: Looking at a Single BGC Float
% In this exercise, we download the NetCDF file for a BGC float deployed
% near Bermuda, inspect its contents, show the trajectory, plot profiles
% for unadjusted and adjusted data, and show the effect of adjustments
% made to the pH concentrations.

% Download NetCDF file for float #5906435
WMO = 5906435;
success = download_float(WMO);
if ~success
    warning('Sprof file for float 5906435 could not be downloaded')
end

% List the sensors for 5906435
list_sensors(WMO);


%% Exercise 2.2: Load the float data
% We see that PH_IN_SITU_TOTAL is available, so load it (along with TEMP and PSAL) from the NetCDF
[data,mdata] = load_float_data(WMO,... % specify WMO number
    {'PSAL','TEMP','PH_IN_SITU_TOTAL'}); % specify variables
fprintf('\nThese fields have been loaded into "data.F%d":\n', WMO)
disp(data.(['F' num2str(WMO)])) % show data that have been loaded into MATLAB
fprintf('\nThese fields have been loaded into "mdata.F%d":\n', WMO)
disp(mdata.(['F' num2str(WMO)])) % show meta data that have been loaded into MATLAB


%% Exercise 2.3: Show the trajectory of the downloaded float with estimated values
show_trajectories(WMO, 'mark_estim', 'yes', 'title', ...
    'Float Trajectory');


%% Exercise 2.4: Show all profiles for salinity and pH from the downloaded float
% this plots the raw, unadjusted data, and includes multiple profiles.
show_profiles(WMO, {'PSAL';'PH_IN_SITU_TOTAL'},'obs','on','raw','yes');


%% Exercise 2.5: Same as above but omitting bad (qc flag 4) data
show_profiles(WMO, {'PH_IN_SITU_TOTAL'},'obs','on','raw','yes','qc',[1 2 3],...
    'title_add','QC Flag 4 Removed');

% Show adjusted profiles with good (flag 1) and probably-good (flag 2) data
show_profiles(WMO, {'PH_IN_SITU_TOTAL'},'obs','on','qc',[1 2],...
    'title_add','Adjusted, QF: 1 and 2');


%% Exercise 2.6: Show sections of raw pH
% This shows the raw, unadjusted data (pcolor plot)
% Mixed layer depth is shown based on the temperature threshold
% (set the value to 2 after 'mld' to use the density threshold instead)
show_sections(WMO, {'PH_IN_SITU_TOTAL'},...
    'mld', 1,...  % tells the function to plot mixed layer depth using T
    'raw','yes'); % tells the function to plot raw data


%% Exercise 2.7: Show sections of adjusted pH
show_sections(WMO, {'PH_IN_SITU_TOTAL'},...
    'mld', 2,...  % tells the function to plot mixed layer depth using rho
    'raw','no',... % tells the function to plot adjusted data if available
    'qc',[1 2]); % tells the function to plot good and probably-good data


%% Exercise 2.8: Show sections of adjusted pH with data restriction
% Since there are no good data after Oct 2024, restrict the plot in time
show_sections(WMO, {'PH_IN_SITU_TOTAL'}, 'mld', 2, 'end', [2024,10,01], ...
    'time_label', 'y', 'qc',[1 2]); % plot only good and probably-good data


%% End of Exercise 2: Clean up the workspace
fprintf(['\nYou have finished Exercise 2! \nDuring the breakout portion, you can ',...
    'go back to this exercise and practice \nwith different WMOs, plot different ',...
    'variables, and explore the additional \narguments and larger capabilities ',...
    'of the functions demonstrated.\n\nProceeding will clear some variables and ',...
    'close any open plots.\n\n'])

do_pause();

clear data mdata success WMO ans
clc; close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 3.1: Oxygen from BGC-Floats around HOT
% In this exercise, we define a region near Hawaii along with a duration of
% time. We identify the float profiles matching those criteria, show
% their trajectories, plot all the matching profiles on one figure, and
% show sections for the unadjusted and adjusted values of dissolved oxygen.


%% Exercise 3.2: Set lat, lon and time limits near Hawaii
latlim=[22 26];
lonlim=[-160 -155];
t1=[2000 1 1];
t2=[];


%% Exercise 3.3: Select floats and profiles based on those limits
[HW_floats,HW_float_profs] = select_profiles(lonlim,latlim,t1,t2,...
    'outside','none', ... % exclude profiles outside the time and space limits
    'type','bgc', ... % include only BGC floats
    'dac','aoml'); % include only floats maintained at AOML

% display the number of matching floats and profiles
disp(' ');
disp(['# of matching profiles: ' num2str(sum(cellfun('length',...
    HW_float_profs)))]);
disp(['# of matching floats: ' num2str(length(HW_floats))]);
disp(' ');


%% Exercise 3.4: Show trajectories for the matching floats, along with the geo limits
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


%% Exercise 3.5: Show trajectories for the matching profiles from each float, along with the geo limits
% Adding the optional input of 'float_profs' with the per-float profile numbers given by
% the select_profiles function will plot only the locations of those
% profiles from the specified floats
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


%% Exercise 3.6: Show sections oxygen for the fifth float in the list of Hawaii floats
% This shows the adjusted data, which is the default setting.
show_sections(HW_floats(5), {'DOXY'}, 'mld', 1);


%% Exercise 3.7: Let's demonstrate adding additional argument by setting a depth range:
% As seen in the plot above, not all floats have data to the full 2000 m
show_sections(HW_floats(5), {'DOXY'}, 'mld', 1, 'depth',[0 1000]);

%% Exercise 3.8 Show good profiles of the section created in Exercise 3.7
% In the previous two plots we saw many missing data. Now we identify the
% oxygen profiles of this float that have good data and show the section
% for them.
% First we load the data into memory so that we can take a closer look
Data_HW5 = load_float_data(HW_floats(5), 'DOXY');
% Then we find all profiles that have at least one good data point
str_floatnum5 = sprintf('F%d', HW_floats(5));
good_doxy5 = find(any(isfinite(Data_HW5.(str_floatnum5).DOXY_ADJUSTED)));
% Plot all profiles with good data
show_sections(HW_floats(5), {'DOXY'}, 'mld', 1, 'depth',[0 1000], ...
    'float_profs', {good_doxy5});

%% Exercise 3.8: Show time series of near-surface oxygen for two Hawaii floats
% Show both floats in one plot per variable, use adjusted values
show_timeseries(HW_floats(4:5), {'DOXY'}, 20, 'per_float', 0);

%% Exercise 3.9: Show good profiles in the timeseries in Exercise 3.8
% So we find the time range of good data for both floats.
Data_HW4 = load_float_data(HW_floats(4), 'DOXY');
str_floatnum4 = sprintf('F%d', HW_floats(4));
good_doxy4 = find(any(isfinite(Data_HW4.(str_floatnum4).DOXY_ADJUSTED)));

% determine first and last date for which either float has good data
min_date = min([Data_HW4.(str_floatnum4).TIME(1,good_doxy4(1)), ...
    Data_HW5.(str_floatnum5).TIME(1,good_doxy5(1))]);
max_date = max([Data_HW4.(str_floatnum4).TIME(1,good_doxy4(end)), ...
    Data_HW5.(str_floatnum5).TIME(1,good_doxy5(end))]);

% show the timeseries with the reduced time range
show_timeseries(HW_floats(4:5), {'DOXY'}, 20, 'per_float', 0, ...
    'start', [year(min_date), month(min_date), day(min_date)], ...
    'end', [year(max_date), month(max_date), day(max_date)]);


%% End of Exercise 3: Clean up the workspace
fprintf(['\nYou have finished Exercise 3! \nDuring the breakout portion, you can ',...
    'go back to this exercise and practice \nwith different WMOs, different ',...
    'lat/lon and time contraints, plot different \n',...
    'variables, and explore the additional arguments and larger capabilities ',...
    'of \nthe functions demonstrated.\n\nProceeding will clear some variables and ',...
    'close any open plots.\n\n'])

do_pause()

clear lonlim latlim t1 t2 ans HW_floats HW_float_profs Data_HW4 Data_HW5
clc; close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 4.1: Nitrate from BGC Floats around BATS
% In this exercise, we define a region near BATS along with a duration of
% time. We identify the float profiles matching those criteria, show
% plot the profiles of each float, look at the time series at varying
% depths with each float, and investigate interesting floats from the time
% series using a map display of the measurements at depth.

%% Exercise 4.2: Set lat, lon and time limits near Bermuda
latlim=[26 36];
lonlim=[-70 -60];
t1=[2000 1 1];
t2=[];


%% Exercise 4.3: Select floats and profiles based on those limits
[BR_floats,BR_float_profs] = select_profiles(lonlim,latlim,t1,t2,...
    'outside','none', ... % exclude profiles outside the time and space limits
    'type','bgc', ... % include only BGC floats
    'dac','aoml',... % include only floats maintained at AOML
    'mode','AD',... % QC mode: A - Adjusted, D - Delayed Mode Quality Control
    'sensor','NITRATE'); % floats must have Nitrate sensor


% display the number of matching floats and profiles
disp(' ');
disp(['# of matching profiles: ' num2str(sum(cellfun('length',...
    BR_float_profs)))]);
disp(['# of matching floats: ' num2str(length(BR_floats))]);
disp(' ');


%% Exercise 4.4: Show profiles of nitrate from all the floats
% Similar to show profiles in the first exercise, but this time with a
% group of floats.
% 'per_float' controls if you make a plot per float (1) or a plot with all
% profiles printed (0).
show_profiles(BR_floats, {'NITRATE'},'per_float',0,'obs','on','qc',[1 2],...
    'float_profs', BR_float_profs, 'title_add','; Adjusted, QF: 1 and 2');


%% Exercise 4.5: Show time series of nitrate at different depths
% Show floats in one plot per variable, use adjusted values
show_timeseries(BR_floats, {'NITRATE'}, [50, 500], ...
    'per_float', 0, 'qc',[1 2]);


%% Exercise 4.6: Let's took a look at a float that stood out in the timeseries
WMO = 5906438; % WMO was pulled from the legend in the last plot
show_trajectories(WMO, 'mark_estim', 'yes', 'title', ...
    'Float trajectory');


%% Exercise 4.7: Nitrate Map
% Show nitrate maps at 50 and 500 dbar depths including all selected floats
show_maps(BR_floats,{'NITRATE'},[50,500],...
    'float_profs',BR_float_profs,'qc',[1 2]);

% Show nitrate maps at 50 and 500 dbar depths for the float of interest
% defined in Exercise 4.6
show_maps(WMO,{'NITRATE'},[50,500],'qc',[1 2]);


% ******** Nitrate around BATS - More Data Visualization Tools Available!
% The Nitrate story around BATS is super interesting, and floats can help
% dig into that cool story.  Stay tuned for another demo on a different
% bgc-argo float visualization tool, webODV for BGC-Argo floats!


%% End of Exercise 4: Clean up the workspace
fprintf(['\nYou have finished Exercise 4! \nDuring the breakout portion, you can ',...
    'go back to this exercise and practice \nwith different WMOs, different ',...
    'lat/lon and time contraints, plot different \n',...
    'variables, and explore the additional arguments and larger capabilities ',...
    'of \nthe functions demonstrated.\n\nProceeding will clear some variables and ',...
    'close any open plots.\n\n'])

do_pause()

clear WMO latlim lonlim t1 t2 ans BR_float_profs BR_floats
clc; close all



%%%%%%%%%%%%%%%%%%%%%%%%%%% TOOLBOX BREAKOUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Suggestions for WMO's, and different lat/lon/time ranges to try!

% Want to look at some more singular floats?
% Interesting floats to look at:
% 5901468 -> first MBARI float with NO3 sensor, deployed by HOT
% 5904859 -> SOCCOM float with many under-ice seasons (Check out it's
%           trajectory!)

% Website of SOCCOM floats - lists SOCCOM WMOs:
% https://www3.mbari.org/soccom/tables/SOCCOM_float_performance.html

% Website of GOBGC floats - lsit GOBGC WMOs:
% https://www3.mbari.org/gobgc/tables/GOBGC_float_performance.html

% Want to look at some floats over a different spatial range?

% Southern Ocean Pacifc Basin ~Lat and ~Lon bounds:
% latlim=[-30 -75];
% lonlim=[115 -70];
% t1=[2014 1 1]; % Dates from 2014 on include lots of data, thanks SOCCOM!
% t2=[];

% Southern Ocean Atlantic Basin Lat and Lon bounds:
% latlim=[-30 -75];
% lonlim=[-70 20];
% t1=[2014 1 1]; % Dates from 2014 on include lots of data, thanks SOCCOM!
% t2=[];

% Southern Ocean Indian Basin Lat and Lon bounds:
% latlim=[-30 -75];
% lonlim=[20 115];
% t1=[2014 1 1]; % Dates from 2014 on include lots of data, thanks SOCCOM!
% t2=[];

% Drake Passage Lat and Lon bounds:
% latlim=[-52 -63];
% lonlim=[-70 -50];
% t1=[2014 1 1]; % Dates from 2014 on include lots of data, thanks SOCCOM!
% t2=[];

% Indian Ocean Lat and Lon bounds:
% latlim=[15 -30];
% lonlim=[20 115];
% t1=[2022 10 14]; % First profile of the first full BGC Indian Ocean float was Oct. 15th 2022!
% t2=[];

% There are many options for the select_profiles options.
% Use "help select_profiles" to list them all (under OPTIONAL INPUTS).
% Add these options after the "t2" argument.
% [floats, float_profs] = select_profiles(lonlim,lat_lim,t1,t2);


%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTRA EXERCISES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Extra Exercise 1.1: Looking at the full netcdf file of a float
% If you wanted to look directly at the netcdf being used rather than using
% a function like list_sensors to examine the netcdf, you are able still
% able to.

% Download NetCDF file for float #5906435
WMO = 5906435;
success = download_float(WMO);
if ~success
    warning('Sprof file for float 5906435 could not be downloaded')
end


%% Extra Exercise 1.2: Display attributes, dimensions, and variables available in the NetCDF
fprintf('\nA look at the Sprof netcdf file of float %d:\n', WMO)
ncdisp(['./Profiles/' num2str(WMO) '_Sprof.nc']) % show the full nc file


%% Extra Exercise 1.3: Extract informational data from the NetCDF
S = ncinfo(['./Profiles/' num2str(WMO) '_Sprof.nc']);
fprintf('\nThese variables are available for float %d:\n', WMO)
{S.Variables.Name}' % show the available variables
% The list_sensors function is designed to do this, just in fewer and
% cleaner steps!



%%%%%%%%%%%%%%%%%%%%%%%%%%%% MORE TUTORIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Other Tutorial Style Scripts Available to You:
% Downloading the OneArgo Toolbox gave you access to a few other tutorials!

% If you want to step through the other existing tutorials, they are:
% Tutorial.m % 2021 GO-BGC Workshop Tutorial
% Show_deep_floats.m % select deep floats with oxygen sensor
% Show_oxygen_maps.m % plots maps for oxygen data that have been DMQC-ed
% Show_under_ice.m % shows how to plot estimated positions under ice


