%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show_under_ice.m
%
% This example script for the OneArgo-Mat toolbox shows the handling
% of missing and estimated positions, which often (but not only) occur
% when a float is under ice.
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

%%
% Floats cannot surface when they are under sea ice.
% Profile data are logged and transmitted when the float can surface again.
% In many cases, the Sprof files from the GDAC will have estimated
% positions for the profiles that were taken under sea ice.
% These are marked by a value of 8 in the POSITION_QC variable.

float1 = 5904859; % this float has under ice interpolation in the Sprof file

% this is the default way of plotting a trajectory - all locations
% are shown in the same color:
show_trajectories(float1, 'title', 'Trajectory partially under ice: Example 1');

% optionally, the estimated locations can be shown in light gray
% (or a different color, based on the value used for
% Settings.color_estim_loc, which is set in initialize_argo):
show_trajectories(float1, 'title', 'Trajectory partially under ice: Example 1', ...
    'mark_estim', 'yes');

%%
% For some floats, under ice locations are marked as missing with a value
% of 9 in the POSITION_QC. By default, the missing locations are no longer
% linearly interpolated between good location values. 
% (Note that this is a change from earlier version of this toolbox.
% Now the default is set in initialize_argo.m with the assignment to
% Settings.interp_lonlat ('no' by default).
% A warning message with the number of interpolated positions
% (13 in this case) is printed in the command window.
% Note that the Sprof file is not modified at all - if the interpolation
% of values should not be used, this can be done by adding optional
% arguments 'interp_lonlat','no' to the calls to show_trajectories and
% load_float_data.

float2 = 4901040;

% the interpolation is by default suppressed - note that in this case no
% locations are shown in gray, only the first and last location are shown
% since these are the only ones with known good values:
show_trajectories(float2, 'title', 'Trajectory partially under ice: Example 2');

% missing locations can be estimated by passing
% 'interp_lonlat', 'yes'
% to the show_trajectories call
% all locations are shown in the same color again in this case:
show_trajectories(float2, 'title', 'Trajectory partially under ice: Example 2', ...
    'interp_lonlat', 'yes');

% again, the estimated locations can be shown in gray - this is done
% in exactly the same way as for locations that were estimated in the
% Sprof file already:
show_trajectories(float2, 'title', 'Trajectory partially under ice: Example 2', ...
    'mark_estim', 'yes', 'interp_lonlat', 'yes');

