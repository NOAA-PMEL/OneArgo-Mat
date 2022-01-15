%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show_under_ice.m
%
% This example script for the BGC-Argo MATLAB toolbox shows the handling
% of missing and estimated positions, which often (but not only) occur
% when a float is under ice.
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
% of 9 in the POSITION_QC. By default, the missing locations are 
% linearly interpolated between good location values. A warning message
% with the number of interpolated positions (13 in this case) is printed
% in the command window.
% Note that the Sprof file is not modified at all - if the interpolation
% of values should not be used, this can be done by adding optional
% arguments 'interp_lonlat','no' to the calls to show_trajectories and
% load_float_data.

float2 = 4901040;

% in the default way of plotting a trajectory all locations
% are shown in the same color again, as in the previous case:
show_trajectories(float2, 'title', 'Trajectory partially under ice: Example 2');

% again, the estimated locations can be shown in gray - this is done
% in exactly the same way as for locations that were estimated in the
% Sprof file already:
show_trajectories(float2, 'title', 'Trajectory partially under ice: Example 2', ...
    'mark_estim', 'yes');

% the interpolation can be suppressed - note that in this case no
% locations are shown in gray, only the first and last location since
% these are the only ones with known good values:
show_trajectories(float2, 'title', 'Trajectory partially under ice: Example 2', ...
    'mark_estim', 'yes', 'interp_lonlat', 'no');

