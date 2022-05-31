function initialize_prof(file_name)
% initialize_prof  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   initialize_prof(file_name)
%
% DESCRIPTION:
%   This function initializes the global struct Prof by reading
%   the index file and processing its information.
%
% INPUTS:
%   file_name : name of the index file (with local path)
%
% OUTPUTS: None. Prof is filled with fields.
%
% AUTHORS:
%   H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL), N. Buzby (UW)
%
% CITATION:
%   H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2022. OneArgo-Mat:
%   A MATLAB toolbox for accessing and visualizing Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.6588042
%
% LICENSE: oneargo_mat_license.m
%
% DATE: JUNE 1, 2022  (Version 1.0.1)

global Prof;

fid = fopen(file_name);
H = textscan(fid,'%s %s %f %f %s %d %s %s','headerlines',9,...
    'delimiter',',','whitespace','');
fclose(fid);
Prof.urls = H{1};
Prof.date = H{2};
Prof.lat  = H{3};
Prof.lon  = H{4};
Prof.ocean = H{5};
Prof.profiler = H{6}; % profiler type
% column 7: institution
Prof.update = H{8};

% the split_sens field is needed by select_profiles_per_type
pT = {'PRES';'TEMP'}; % for old floats without salinity sensor
pTS = {'PRES';'TEMP';'PSAL'}; % for all new core and deep floats
Prof.split_sens = cell(length(Prof.urls), 1);
Prof.split_sens(:) = {pTS};
Prof.split_sens(Prof.profiler == 845) = {pT};

% adjust longitude to standard range of -180..180 degrees
Prof.lon(Prof.lon > 180) = Prof.lon(Prof.lon > 180) - 360;
Prof.lon(Prof.lon < -180) = Prof.lon(Prof.lon < -180) + 360;

Prof.wmo = extractBetween(Prof.urls, '/', '/');
