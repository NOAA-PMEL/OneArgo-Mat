function datenum_value = check_datenum(date_value)
% check_datenum  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   datenum_value = check_datenum(date_value)
%
% DESCRIPTION:
%   This function checks if the given date_value conforms to the
%   datenum format. If so, the corresponding datenum_value is returned.
%   If not, a warning is issued and an empty array is returned.
%
% INPUT:
%   date_value : date (array) in one of the following formats:
%                [YYYY MM DD HH MM SS] or [YYYY MM DD]
%
% OUTPUT:
%   datenum_value : the corresponding datenum value (day 1 is Jan-1-0000),
%                e.g., datenum([2020,1,1]) is 737791
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

if isnumeric(date_value) && ...
        (length(date_value) == 3 || length(date_value) == 6)
    datenum_value = datenum(date_value);
else
    warning(['dates should be in [YYYY MM DD HH MM SS] or ', ...
        '[YYYY MM DD] format']);
    datenum_value = [];
end
