function [main_sensor, sensor_number] = get_sensor_number(sensor_name)
% get_sensor_number  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   [main_sensor, sensor_number] = get_sensor_number(sensor_name)
%
% DESCRIPTION:
%   This function returns the name of the main sensor and, if applicable,
%   the number of the sensor. If sensor_name represents the main
%   sensor, sensor_number is an empty string, not '1'. sensor_number
%   is only defined for additional sensors, e.g. '2' for 'DOXY2'.
%
% INPUTS:
%   sensor_name : the full name of the sensor, e.g.: 'TEMP' or 'DOXY2'
%
% OUTPUTS:
%   main_sensor : the name of the main sensor, e.g. 'TEMP' or 'DOXY'
%   sensor_number : the number of the sensor returned as a string,
%                 e.g.: '2' (empty string for the main sensor)
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

% BSDOXY FIX: pristine v1.1.0 could not parse the underscore form 'DOXY_2'
% (it expects 'DOXY2') and crashed on unknown sensors. Parse the base sensor
% name + optional sensor number, handling DOXY, DOXY2, DOXY_2, BBP700, BBP700_2,
% DOWNWELLING_PAR, PH_IN_SITU_TOTAL, ... and return empty for a genuinely
% unrecognised name (initialize_sprof expects empty to flag it).
%   main : [A-Z]+ , optional _LETTER groups (DOWNWELLING_PAR) , optional 3-digit
%          wavelength (BBP700).   num : optional trailing 2..9, with or without '_'.
match = regexp(sensor_name, ...
    '^(?<main>[A-Z]+(_[A-Z]+)*(\d{3})?)(?<num>_?[2-9])?$', 'names');
if isempty(match)
    main_sensor = '';
    sensor_number = '';
    return
end
main_sensor = match.main;
sensor_number = strrep(match.num, '_', '');
