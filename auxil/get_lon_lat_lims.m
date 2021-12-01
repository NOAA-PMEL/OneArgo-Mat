function [lon_lim, lat_lim, Data] = get_lon_lat_lims(Data)
% get_lon_lat_lims  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   [lon_lim, lat_lim, Data] = get_lon_lat_lims(Data)
%
% DESCRIPTION:
%   This function obtains maximum and minimum latitude and longitude values
%   from input data.
%
% PREREQUISITE: 
%   Sprof file(s) for the specified float(s) must exist locally.
%
% INPUT:
%   Data     : a struct that must contain at least LONGITUDE and LATITUDE
%              fields
%
% OUTPUTS:
%   lon_lim  : a 2-element vector with minimum,maximum longitudes
%              (normally using a range of -180..180 degrees)
%   lat_lim  : a 2-element vector with minimum,maximum latitudes
%   Data     : if using a range of 0..360 degrees instead of the default
%              range of -180..180 degrees results in a shorter range to
%              cover all trajectories, a field ALT_LON is added, which uses
%              a range of 0..360 degrees;
%              it is unchanged from the input otherwise
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

floats = fieldnames(Data);
nfloats = length(floats);

% a range of -180..180 degrees is the standard convention for longitude
lon_lim = [180, -180];
lat_lim = [90, -90];
% the alternate range of 0..360 degrees is tested as well
lon_lim2 = [360, 0];

for i = 1:nfloats
    lon_lim(1) = min([lon_lim(1), min(Data.(floats{i}).LONGITUDE(:))]);
    lon_lim(2) = max([lon_lim(2), max(Data.(floats{i}).LONGITUDE(:))]);
    lat_lim(1) = min([lat_lim(1), min(Data.(floats{i}).LATITUDE(:))]);
    lat_lim(2) = max([lat_lim(2), max(Data.(floats{i}).LATITUDE(:))]);
    ALT_LON{i} = Data.(floats{i}).LONGITUDE;
    ALT_LON{i}(ALT_LON{i} < 0) = ALT_LON{i}(ALT_LON{i} < 0) + 360;
    lon_lim2(1) = min([lon_lim2(1), min(ALT_LON{i})]);
    lon_lim2(2) = max([lon_lim2(2), max(ALT_LON{i})]);
end

if diff(lon_lim2) < diff(lon_lim)
    % the alternate (0..360 degrees) longitude range is shorter; use it
    lon_lim = lon_lim2;
    for i = 1:nfloats
        Data.(floats{i}).ALT_LON = ALT_LON{i};
    end
end
