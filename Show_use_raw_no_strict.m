%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show_use_raw_no_strict.m
%
% In this example script for the OneArgo-Mat toolbox
% we define a region in the Northeast Pacific along with
% a duration of time, and identify the float profiles matching those
% criteria. We show the trajectories of all the matching floats and plot
% profiles that match the criteria for one of the floats.
% It shows the use of the 'raw','no_strict' option for plotting
% sections, profiles, and time series. Only profiles in adjusted
% mode will be shown if this option is used.
%
% AUTHORS:
%   H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL), N. Buzby (UW)
%
% CITATION:OSP
%   H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2022. OneArgo-Mat:
%   A MATLAB toolbox for accessing and visualizing Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.6588041
%
% LICENSE: oneargo_mat_license.m
%
% DATE: JUNE 1, 2022  (Version 1.0.1)

%% Close figures, clean up workspace, clear command window
close all; clear; clc

%% Initialize
% This function defines standard settings and paths and creates Index
% and Profiles folders in your current path if necessary. 
% It also downloads the relevant Sprof index file from the GDAC to your 
% Index folder. The Sprof index is referenced when downloading and
% subsetting float data based on user specified criteria in other functions.
initialize_argo();

%% Set limits in the NE Pacific from 2016 to 2020
latlim=[45 65];
lonlim=[-170 -115];
t1=[2010 1 1];
t2=[2020 12 31];

%% select profiles based on those limits with specified sensor (DOXY)
[NEP_floats,NEP_float_profs] = select_profiles(lonlim,latlim,t1,t2,...
    'sensor','DOXY',... % this selects only floats with oxygen sensors
    'mode','D');

% display the number of matching floats and profiles
disp(' ');
disp(['# of matching profiles: ', ...
    num2str(sum(cellfun('length', NEP_float_profs)))]);
disp(['# of matching floats: ' num2str(length(NEP_floats))]);
disp(' ');

%% Show trajectories for the matching floats
% This function loads the data for plotting.
% Adding the optional input pair 'color','multiple' will plot different
% floats in different colors.
show_trajectories(NEP_floats,'color','multiple');

%%
% In this setup, the third float (5900961) doesn't have DOXY_ADJUSTED
% values, so - to show the same variable type (raw vs adjusted) for
% all floats, raw valuew are shown.
% However, raw DOXY values should normally not have QC flags of 1 or 2.
% (Float 4901137 still does at this point, but this may/should change.)
% For the other two floats, there are no "good raw" DOXY values,
% so the plot isn't shown.
show_sections(NEP_floats(1:3),{'DOXY'}, ...
    'float_profs',NEP_float_profs(1:3),...
    'raw','no',...
    'qc',[1 2]);
%%
% In this setup, the plot is shown for DOXY_ADJUSTED for the first
% two floats, which have adjusted oxygen values.
show_sections(NEP_floats(1:3),{'DOXY'}, ...
    'float_profs',NEP_float_profs(1:3),...
    'raw','no_strict',...
    'qc',[1 2]);

%%
% Switches to 'raw' mode for all floats. With the given QC flags,
% plots are empty for at least two floats.
show_profiles(NEP_floats(1:3),{'DOXY'}, ...
    'float_profs',NEP_float_profs(1:3),...
    'per_float',1,...
    'raw','no',...
    'qc',[1 2]);

%%
% Again, plots will only be shown for the first two floats.
show_profiles(NEP_floats(1:3),{'DOXY'}, ...
    'float_profs',NEP_float_profs(1:3),...
    'per_float',1,...
    'raw','no_strict',...
    'qc',[1 2]);

%%
% Switches to raw mode for all floats.
show_timeseries(NEP_floats(1:3),{'TEMP'}, ...
    990,...
    'float_profs',NEP_float_profs(1:3),...
    'per_float',1,...
    'raw','no',...
    'var2','DOXY',...
    'qc',[1 2]);

%%
% Shows plots for the two floats with DOXY_ADJSUTED values.
show_timeseries(NEP_floats(1:3),{'TEMP'}, ...
    990,...
    'float_profs',NEP_float_profs(1:3),...
    'per_float',1,...
    'raw','no_strict',...
    'var2','DOXY',...
    'qc',[1 2]);
