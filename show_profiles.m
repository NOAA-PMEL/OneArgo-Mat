function good_float_ids = show_profiles(profile_ids, variables, varargin)
% show_profiles  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   good_float_ids = show_profiles(profile_ids, variables, varargin)
%
% DESCRIPTION:
%   This an intermediary function that downloads profile(s) for the given
%   float(s) and calls plot_profile to create the plot(s).
%
% INPUTS:
%   profile_ids  : internally used indices of individual profiles
%   variables    : cell array of variable(s) (i.e., sensor(s)) to show 
%                  (if not set: {'DOXY'} (=O2) is used)
%
% OPTIONAL INPUTS:
%  'type',type   : by default (type='profiles'), the given IDs refer to
%                  profile IDs (obtained with select_profiles); use
%                  'type','floats' to show the profiles of a given float
%  'obs',on/off  : by default (obs='off') only lines are shown for each
%                  profile; 'obs','on' shows points on the profile at which
%                  each measurement was made
%  'raw',raw     : plot raw, i.e., unadjusted data if set to 'yes';
%                  default: 'no' (i.e., plot adjusted data if available)
%  'qc',flags    : show only values with the given QC flags (as an array)
%                  0: no QC was performed; 
%                  1: good data; 
%                  2: probably good data;
%                  3: probably bad data that are potentially correctable;
%                  4: bad data; 
%                  5: value changed; 
%                  6,7: not used;
%                  8: estimated value; 
%                  9: missing value
%                  default setting: 
%                  [1,2] for adjusted data; [0:9] for raw data
%                  See Table 7 in Bittig et al.:
%                  https://www.frontiersin.org/files/Articles/460352/fmars-06-00502-HTML-r1/image_m/fmars-06-00502-t007.jpg
%  'title_add',text : add the given text to all titles
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

if nargin < 2
    warning('Usage: show_profiles(profile_ids, variables, varargin)')
end

if isempty(profile_ids)
    warning('no profiles specified')
    return
end

% set defaults
if nargin < 2
    variables = {'DOXY'};
end
if ~nargin
    profile_ids = Settings.demo_float;
end
type = 'profiles';
varargpass= {};

% parse optional arguments
for i = 1:2:length(varargin)
    if strcmpi(varargin{i}, 'type')
        type = varargin{i+1};
    else
        varargpass = [varargpass, varargin{i:i+1}];
        if strcmpi(varargin{i}, 'qc')
            if min(varargin{i+1}) < 0 || max(varargin{i+1}) > 9
                warning('only QC flags 0..9 are allowed!')
            end
        end
    end
end

% convert requested variable to cell array if necessary (string was used)
if ischar(variables)
    variables = cellstr(variables);
end

if strncmpi(type, 'prof', 4)
    % profile IDs need to be converted to float IDs
    all_float_ids = str2num(cell2mat(Sprof.wmo(profile_ids)));
else
    all_float_ids = profile_ids;
end
uniq_float_ids = unique(all_float_ids);

% download Sprof files if necessary
good_float_ids = download_multi_floats(uniq_float_ids);

if isempty(good_float_ids)
    warning('no valid floats found')
else
    if strncmpi(type, 'prof', 4)
        for i = 1:length(good_float_ids)
            idx = (all_float_ids == good_float_ids(i));
            float_profs{i} = Sprof.fprofid(profile_ids(idx));
        end
        [Data, Mdata] = load_float_data(good_float_ids, variables, ...
            float_profs);
    else
        [Data, Mdata] = load_float_data(good_float_ids, variables);
    end
    plot_profiles(Data, Mdata, variables, varargpass{:});
end

