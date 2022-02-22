%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show_missing_lonlat.m
%
% This example script for the BGC-Argo MATLAB toolbox shows the revised
% handling of missing positions, which often (but not only) occur
% when a float is under ice.
% This is a change in behavior from versions 1.0 and 1.1 of the toolbox.
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
% DATE: FEBRUARY 22, 2022  (Version 1.2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Illustration of behavior for versions 1.0 and 1.1 of the toolbox
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

%% Illustration of behavior for version 1.2 of the toolbox
% by default all positions (known and estimated) are shown in the same
% color
show_trajectories(float, 'title', ...
    'All positions of F4901040');

% if the 'mark_estim','yes' option is used, estimated positions are
% are shown in gray
show_trajectories(float, 'mark_estim','yes', 'title', ...
    'Known and estimated positions of F4901040');

% if interpolation is not explicitly suppressed (see above), profiles
% with estimated locations are included in the search
[floats,float_profs] = select_profiles([175,185],[80,85],...
    [2008,10,1],[2009,1,1])
