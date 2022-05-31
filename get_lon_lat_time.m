function [lon, lat, time] = get_lon_lat_time(float_ids, float_profs)
% get_lon_lat_time  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   [lon, lat, time] = get_lon_lat_time(float_ids [, float_profs])
%
% DESCRIPTION:
%   This function loads longitude, latitude, and time information
%   for the specified floats (and their specified profiles, if given).
%
% INPUT:
%   float_ids   : WMO ID(s) of one or more floats
%
% OPTIONAL INPUT:
%   float_profs : cell array with indices of selected profiles (per float,
%                 not global)
%
% OUTPUTS:
%   lon  : cell array with longitude values for all specified floats
%   lat  : cell array with latitude values for all specified floats
%   time : cell array with time values for all specified floats (in
%          datenum format)
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

if isempty(float_ids)
    warning('no floats specified');
    lon = [];
    lat = [];
    time = [];
    return
end

if nargin < 2
    float_profs = [];
end

Data = load_float_data(float_ids, {}, float_profs);
nfloats = length(fieldnames(Data));
lon = cell(nfloats, 1);
lat = cell(nfloats, 1);
time = cell(nfloats, 1);

for f = 1:nfloats
    str_floatnum = ['F', num2str(float_ids(f))];
    lon{f} = Data.(str_floatnum).LONGITUDE(1,:)';
    lat{f} = Data.(str_floatnum).LATITUDE(1,:)';
    time{f} = Data.(str_floatnum).TIME(1,:)';
end
