function good_float_ids = show_trajectories(float_ids,varargin)
% show_trajectories  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   show_trajectories(float_ids,varargin)
%
% DESCRIPTION:
%  This is an intermediary function that downloads profiles for at least
%   one given float and calls plot_trajectories to create the plot.
%
% INPUTS:
%   float_ids : WMO ID(s) of one or more floats 
%               (if not set: Settings.demo_float is used as a demo)
%
% OPTIONAL INPUTS:
%   'color',color : color (string) can be either 'multiple' (different
%                   colors for different floats), or any standard Matlab
%                   color descriptor ('r', 'k', 'b', 'g' etc.)
%                   (all trajectories will be plotted in the same color)
%                   default value is 'r' (red)
%  'float_profs',fp : fp is an array with the per-float indices of the
%                   selected profiles, as returned by function
%                   select_profiles - use this optional argument if you
%                   don't want to plot the full trajectories of the
%                   given floats, but only those locations that match
%                   spatial and/or temporal constraints
%  'png',fn_png   : save the plot to a png file with the given 
%                   file name (fn_png)
%  'title',title  : title for the plot (default: "Float trajectories")
%
% OUTPUT:
%   good_float_ids : array of the float IDs whose Sprof files were
%                    successfully downloaded or existed already
%
% AUTHORS: 
%   H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL),
%   J. Plant, T. Maurer, Y. Takeshita (MBARI), D. Nicholson (WHOI),
%   and A. Gray (UW)
%
% CITATION:
%   H. Frenzel*, J. Sharp*, A. Fassbender, J. Plant, T. Maurer,
%   Y. Takeshita, D. Nicholson, A. Gray, 2021. BGC-Argo-Mat: A MATLAB
%   toolbox for accessing and visualizing Biogeochemical Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.4971318.
%   (*These authors contributed equally to the code.)
%
% LICENSE: bgc_argo_mat_license.m
%
% DATE: June 15, 2021

global Settings;

% make sure Settings is initialized
if isempty(Settings)
    initialize_argo();
end

if ~nargin
    float_ids = Settings.demo_float;
end

% set defaults
color = 'r'; % red
float_profs = [];
fn_png = [];
title1 = 'Float trajectories';

% parse optional arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'color')
        color = varargin{i+1};
    elseif strcmpi(varargin{i}, 'float_profs')
        float_profs = varargin{i+1};
    elseif strcmpi(varargin{i}, 'png')
        fn_png = varargin{i+1};
    elseif strcmpi(varargin{i}, 'title')
        title1 = varargin{i+1};
    end
end

% download Sprof files if necessary
good_float_ids = download_multi_floats(float_ids);

if isempty(good_float_ids)
    warning('no valid floats found')
else
    % meta data return values and observations are not needed here
    Data = load_float_data(good_float_ids,{},float_profs);   
    plot_trajectories(Data, color, title1, fn_png);
end
