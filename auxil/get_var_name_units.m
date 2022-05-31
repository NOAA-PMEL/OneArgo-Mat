function [long_name, units] = get_var_name_units(short_name)
% get_var_name_units  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   [long_name, units] = get_var_name_units(short_name)
%
% DESCRIPTION:
%   This function returns the long name and the units for the variable
%   with the given short name.
%
% INPUT:
%   short_name : case-sensitive name of a variable as it appears in
%                the Sprof index file, e.g., TEMP or DOXY
%
% OUTPUTS:
%   long_name  : long name of the variable
%   units      : units of the variable
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

if nargin < 1
    warning('Usage: get_var_name_units(short_name)')
    return
end

if contains(short_name,'TEMP')
    long_name = 'Temperature';
    units = '(^oC)';
elseif contains(short_name,'PSAL')
    long_name = 'Salinity';
    units = '(PSU)';
elseif contains(short_name,'DENS')
    long_name = 'Density';
    units = '(kg m^{-3})';
elseif contains(short_name,'PRES')
    long_name = 'Pressure';
    units = '(dbar)';
elseif contains(short_name,'DOXY')
    long_name = 'Dissolved Oxygen';
    units = '(\mumol kg^{-1})';
elseif contains(short_name,'NITRATE')
    long_name = 'Nitrate ';
    units = '(\mumol kg^{-1})';
elseif contains(short_name,'AOU')
    long_name = 'Apparent Oxygen Utilization ';
    units = '(\mumol kg^{-1})';
elseif contains(short_name,'PH_IN_SITU')
    long_name = 'pH';
    units = '';
elseif contains(short_name,'CHLA')
    long_name = 'Chlorophyll-a';
    units = '(mg m^{-3})';
elseif contains(short_name,'CDOM')
    long_name = 'Colored Dissolved Organic Matter';
    units = '(ppb)';
elseif contains(short_name,'BBP')
    long_name = 'Backscatter';
    units = '(m^{-1})';
elseif contains(short_name,'DOWN_IRRADIANCE')
    long_name = 'Downwelling irradiance';
    units = '(W m^{-2} nm^{-1})';
elseif contains(short_name,'UP_RADIANCE')
    long_name = 'Upwelling radiance';
    units = '(W m^{-2} nm^{-1})';
elseif contains(short_name,'DOWNWELLING_PAR')
    long_name = 'Downwelling PAR';
    units = '(\mumol Quanta m^{-2} sec^{-1})';
elseif contains(short_name,'BISULFIDE')
    long_name = 'Bisulfide';
    units = '(\mumol kg^{-1})';
elseif contains(short_name,'TURBIDITY')
    long_name = 'Turbidity';
    units = '(ntu)';
elseif strncmp(short_name, 'CP', 2)
    long_name = 'Particle beam attenuation';
    units = '(m^{-1})';
elseif strncmp(short_name, 'CNDC', 2)
    long_name = 'Electrical conductivity';
    units = '(mhos m^{-1})';
else
    warning('unknown variable')
    long_name = [];
    units = [];
end

% add extensions (wavelength and/or sensor number) where applicable
[main_sensor, sensor_number] = get_sensor_number(short_name);

if contains(main_sensor, 'RADIANCE') || strncmp(main_sensor, 'CP', 2) || ...
        strncmp(main_sensor, 'BBP', 3)
    wavelength = regexp(main_sensor, '\d{3}', 'once', 'match');
    if isempty(wavelength)
        if ~isempty(sensor_number)
            long_name = sprintf('%s (sensor %s)', long_name, sensor_number);
        end
    else
        if isempty(sensor_number)
            long_name = sprintf('%s (%s nm)', long_name, wavelength);
        else
            long_name = sprintf('%s (%s nm; sensor %s)', long_name, ...
                wavelength, sensor_number);
        end
    end
elseif ~isempty(sensor_number)
    long_name = sprintf('%s (sensor %s)', long_name, sensor_number);
end
