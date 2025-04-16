function inpoly = get_inpolygon(lon, lat, lon_lim, lat_lim)
% get_inpolygon  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   inpoly = get_inpolygon(lon, lat, lon_lim, lat_lim)
%
% DESCRIPTION:
%   This function determines the lon/lat points that lie within
%   the given lon/lat limits.
%
% INPUTS:
%   lon      : vector of longitude values
%   lat      : vector of latitude values (must have same length as lon)
%   lon_lim  : longitude limits
%   lat_lim  : latitude limits
%            * Latitude and longitude limits can be input as either
%            two element vectors (e.g., [LAT1 LAT2]) for maximum and minimum
%            limits, or as same-sized vectors with at least 3 elements
%            for vertices of a polygon
%            * Longitude can be input in either the -180 to 180 degrees
%            format or 0 to 360 degrees format
%
% OUTPUT:
%   inpoly   : vector of 0s and 1s (same length as lon and lat)
%
% AUTHORS:
%   H. Frenzel and J. Sharp (UW-CICOES), A. Fassbender (NOAA-PMEL), N. Buzby (UW)
%
% CITATION:
%   H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2025. OneArgo-Mat:
%   A MATLAB toolbox for accessing and visualizing Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.6588041
%
% LICENSE: oneargo_mat_license.m
%
% DATE: APRIL 16, 2025  (Version 1.1.0)

% longitude limits can be specified in -180..180 or 0..360 degree formats,
% or any other 360 degree range that encloses all the desired longitude
% values, e.g., [-20 200])
% [10 350] is NOT equivalent to [-10 10]: the former results in a range
% of 340 degrees, the latter in a range of 20 degrees.
if max(lon_lim) > min(lon_lim) + 360
    % (S)prof.lon values use -180..180 range
    lon_lim(lon_lim < -180) = -180;
    lon_lim(lon_lim > 180) = 180;
else
    floor_lon = floor(min(lon_lim));
    if any(lon < floor_lon)
        lon(lon < floor_lon) = lon(lon < floor_lon) + 360;
    else
        ceil_lon = ceil(max(lon_lim));
        lon(lon > ceil_lon) = lon(lon > ceil_lon) - 360;
    end
end

if numel(lon_lim) == 2 && numel(lat_lim) == 2
    lon_vert = [lon_lim(1) lon_lim(2) lon_lim(2) lon_lim(1)];
    lat_vert = [lat_lim(1) lat_lim(1) lat_lim(2) lat_lim(2)];
else
    lon_vert = lon_lim;
    lat_vert = lat_lim;
end

if length(lon_vert) ~= length(lat_vert)
    warning('lon/lat limit vectors must have the same length')
    inpoly = [];
else
    inpoly = inpolygon(lon, lat, lon_vert, lat_vert);
end
