function good_float_ids = download_multi_floats(float_ids)
% download_multi_floats  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   good_float_ids = download_multi_floats(float_ids)
%
% DESCRIPTION:
%   This function downloads *Sprof.nc files for the specified float(s).
%   A message is shown if profiles for any of these floats could not be
%   downloaded.
%
% INPUT:
%   float_ids : numerical array with WMO ID(s) of the float(s)
%
% OUTPUT:
%   good_float_ids : WMO ID(s) of the float(s) whose Sprof files were downloaded
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
    warning('Usage: download_multi_floats(float_ids)')
    return
end

is_good = ones(length(float_ids), 1);
not_found = '';
count = 0;
for i = 1:length(float_ids)
    if ~download_float(float_ids(i))
        is_good(i) = 0;
        not_found = sprintf('%s %d', not_found, float_ids(i));
        count = count + 1;
        % avoid too long lines in command window display
        if count == 10
            not_found = [not_found, newline];
            count = 0;
        end
    end
end
good_float_ids = float_ids(is_good == 1);
if ~isempty(not_found)
    fprintf('Sprof files could not be downloaded for floats:\n%s\n', ...
        not_found);
end
