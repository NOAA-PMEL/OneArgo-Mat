%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show_missing_lonlat.m
%
% This example script for the OneArgo-Mat toolbox shows the revised
% handling of missing positions, which often (but not only) occur
% when a float is under ice.
% This is a change in behavior from versions 1.0 and 1.1 of the BGC-Argo
% version of the toolbox.
% Note that in the latest version of OneArgo-Mat, the default setting
% for interp_lonlat is defined in initialize_argo.m
% It is initially set to 'no' to speed up profile selection etc.,
% but this can be changed by the user.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Illustration of behavior for versions 1.0 and 1.1 of the BGC-Argo toolbox
% this float has 15 profiles, but only two of them have valid lon/lat values
% (both in the index file and the Sprof file)
float = 4901040;
% this is the way the trajectory was displayed in versions 1.0 and 1.1
% of the toolbox
show_trajectories(float, 'interp_lonlat', 'no', 'title', ...
    'Known positions of F4901040');

% if the time period of interest was restricted like this, the float would
% not have been found with select_profiles:
[floats,float_profs] = select_profiles([175,185],[80,85],...
    [2008,10,1],[2009,1,1],'interp_lonlat','no')

%% Illustration of behavior for versions 1.2 and 1.3 of the BGC-Argo toolbox
% by default all positions (known and estimated) are shown in the same
% color
show_trajectories(float, 'interp_lonlat', 'yes', 'title', ...
    'All positions of F4901040');

% if the 'mark_estim','yes' option is used, estimated positions are
% are shown in gray
show_trajectories(float, 'interp_lonlat', 'yes', 'mark_estim', 'yes', ...
    'title', 'Known and estimated positions of F4901040');

% if interpolation is enabled, profiles
% with estimated locations are included in the search
[floats,float_profs] = select_profiles([175,185],[80,85],...
    [2008,10,1],[2009,1,1], 'interp_lonlat', 'yes')
