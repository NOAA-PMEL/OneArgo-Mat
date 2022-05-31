function xl = set_xlim(start_date, end_date)
% set_xlim  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   xl = set_xlim(start_date, end_date)
%
% DESCRIPTION:
%   This function sets the xlimits of the current plot if at least one
%   of the given dates is not empty.
%
% PREREQUISITE:
%   The plot window must exist already. The x axis must use datenum values.
%
% INPUT:
%   start_date : datenum value for the starting point (may be empty)
%   end_date   : datenum value for the ending point (may be empty)
%
% OUTPUT:
%   xl         : the xlim values (in datenum units, [start end])
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

xl = xlim;
if ~isempty(start_date) || ~isempty(end_date)
    if ~isempty(start_date)
        xl(1) = start_date;
    end
    if ~isempty(end_date)
        xl(2) = end_date;
    end
    xlim(xl);
end
