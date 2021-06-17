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
%  'prof_ids',ids : ids is an array with the global indices of the
%                   selected profiles, as returned by function
%                   select_profiles - use this optional argument if you
%                   don't want to plot the full trajectories of the
%                   given floats, but only those that match spatial
%                   and/or temporal constraints
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

global Settings Sprof;

if ~nargin
    float_ids = Settings.demo_float;
end

% set defaults
color = 'r'; % if set, removes profiles outside time/space constraints
prof_ids = {};

% parse optional arguments
for i = 1:2:length(varargin)
    if strcmpi(varargin{i}, 'color')
        color = varargin{i+1};
    elseif strcmpi(varargin{i}, 'prof_ids')
        prof_ids = varargin{i+1};
    end
end

if isempty(prof_ids)
    float_profs = {};
else
    % convert global profile IDs to individual (per float) profile IDs
    float_profs = cell(length(float_ids),1); % pre-allocate
    all_float_ids = str2num(cell2mat(Sprof.wmo(prof_ids))); % get all float ids
    for i = 1:length(float_ids)
        idx = (all_float_ids == float_ids(i)); % index based on float id
        % obtain profiles of that float to plot
        float_profs{i} = Sprof.fprofid(prof_ids(idx));
    end
end

% download Sprof files if necessary
good_float_ids = download_multi_floats(float_ids);

if isempty(good_float_ids)
    warning('no valid floats found')
else
    % meta data return values and observations are not needed here
    Data = load_float_data(good_float_ids,{},float_profs);   
    plot_trajectories(Data,color);
end
