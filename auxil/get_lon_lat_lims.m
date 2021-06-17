function [lon_lim, lat_lim, Data] = get_lon_lat_lims(Data)
% get_lon_lat_lims  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   [lon_lim, lat_lim, Data] = get_lon_lat_lims(Data)
%
% DESCRIPTION:
%   This function obtains maximum and minimum latitude and longitude values
%   from input data
%
% PREREQUISITE: 
%   Sprof file(s) for the specified float(s) must exist locally.
%
% INPUTS:
%   Data     : a struct that must contain at least LONGITUDE and LATITUDE
%              fields
%
% OUTPUTS:
%   lon_lim  : a 2-element vector with minimum,maximum longitudes
%              (normally using a range of -180..180 degrees)
%   lat_lim  : a 2-element vector with minimum,maximum latitudes
%   Data     : if all points are within 30 degrees of 180W/E, a field
%              ALT_LON is added that uses a 0..360 range instead;
%              it is unchanged from the input otherwise
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

floats = fieldnames(Data);
nfloats = length(floats);

% NOTE: assuming -180..180 convention for longitude
lon_lim = [180, -180];
lat_lim = [90, -90];

for i = 1:nfloats
    lon_lim(1) = min([lon_lim(1), min(Data.(floats{i}).LONGITUDE(:))]);
    lon_lim(2) = max([lon_lim(2), max(Data.(floats{i}).LONGITUDE(:))]);
    lat_lim(1) = min([lat_lim(1), min(Data.(floats{i}).LATITUDE(:))]);
    lat_lim(2) = max([lat_lim(2), max(Data.(floats{i}).LATITUDE(:))]);
end

% special treatment for floats that cross the dateline
if lon_lim(1) <= -150 && lon_lim(2) >= 150
    other_pts = 0;
    for i = 1:nfloats
        other_pts = other_pts + ...
            length(find(Data.(floats{i}).LONGITUDE > -150 & ...
            Data.(floats{i}).LONGITUDE < 150));
    end
    % all points are within 30 degrees of the dateline
    if ~other_pts
        lon_lim = [360, 0];
        for i = 1:nfloats
            Data.(floats{i}).ALT_LON = Data.(floats{i}).LONGITUDE;
            Data.(floats{i}).ALT_LON(Data.(floats{i}).ALT_LON < 0) = ...
                Data.(floats{i}).ALT_LON(Data.(floats{i}).ALT_LON < 0) + 360;
            lon_lim(1) = min([lon_lim(1), min(Data.(floats{i}).ALT_LON(:))]);
            lon_lim(2) = max([lon_lim(2), max(Data.(floats{i}).ALT_LON(:))]);
        end
    end
end

