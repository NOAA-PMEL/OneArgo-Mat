function Data = interp_lonlat(Data, floatnum)
% interp_lonlat  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   Data = interp_lonlat(Data, floatnum)
%
% DESCRIPTION:
%   This function interpolates missing location data for one float
%   if it can be done based on good location data.
%
% INPUT:
%   Data     : struct with float data (load_float_data format)
%   floatnum : WMO ID of one float (integer)
%
%
% OUTPUTS:
%   Data     : struct with possible modifications in the
%              LONGITUDE, LATITUDE, POSITION_QC fields for the
%              selected float
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

global Settings;

str_floatnum = ['F', num2str(floatnum)];
% if possible, interpolate lon/lat for profiles whose
% position was marked as missing
% note that profiles with position marked as bad (qc = 4)
% are not affected, but this could be done easily by
% including those profile indices in idx_miss
idx_miss = find(Data.(str_floatnum).POSITION_QC(1,:) == 9);
if idx_miss
    new_values = 0;
    % include qc values of 0,1,2,8 in good values
    idx_good = find(Data.(str_floatnum).POSITION_QC(1,:) < 3 | ...
        Data.(str_floatnum).POSITION_QC(1,:) == 8);
    % need to make sure that no jumps across dateline occur, so
    % use alternate longitude if necessary
    FloatData = struct();
    FloatData.(str_floatnum) = Data.(str_floatnum);
    [~, ~, FloatData] = get_lon_lat_lims(FloatData);
    % use 1D arrays for simplicity
    if isfield(FloatData.(str_floatnum), 'ALT_LON')
        lon = FloatData.(str_floatnum).ALT_LON(1,:);
    else
        lon = Data.(str_floatnum).LONGITUDE(1,:);
    end
    lat = Data.(str_floatnum).LATITUDE(1,:);
    pos_qc = Data.(str_floatnum).POSITION_QC(1,:);
    for i = 2:length(idx_good)
        since_last = idx_good(i-1)+1:idx_good(i)-1;
        idx_miss_here = intersect(since_last, idx_miss);
        if idx_miss_here
            % at least one missing value exists between good values
            lin_intp(idx_good(i-1):idx_good(i)) = ...
                linspace(lon(idx_good(i-1)), ...
                lon(idx_good(i)), ...
                idx_good(i) - idx_good(i-1) + 1);
            lat_intp(idx_good(i-1):idx_good(i)) = ...
                linspace(lat(idx_good(i-1)), ...
                lat(idx_good(i)), ...
                idx_good(i) - idx_good(i-1) + 1);
            lon(idx_miss_here) = lin_intp(idx_miss_here);
            lat(idx_miss_here) = lat_intp(idx_miss_here);
            pos_qc(idx_miss_here) = 8;
            new_values = new_values + length(idx_miss_here);
        end
    end
    if new_values
        if isfield(FloatData.(str_floatnum), 'ALT_LON')
            % revert back to standard range of -180..180 degrees
            lon(lon > 180) = lon(lon > 180) - 360;
        end
        n_levels = size(Data.(str_floatnum).LONGITUDE, 1);
        Data.(str_floatnum).LONGITUDE = repmat(lon, n_levels, 1);
        Data.(str_floatnum).LATITUDE = repmat(lat, n_levels, 1);
        Data.(str_floatnum).POSITION_QC = repmat(pos_qc, n_levels, 1);
        if Settings.verbose
            fprintf(['position interpolation performed for %d ',...
                'missing values of float %d\n'], new_values, floatnum);
        end
    end
end
