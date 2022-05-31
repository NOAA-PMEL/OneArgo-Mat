function good_float_ids = download_tech_files(float_ids)
% download_tech_files  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   good_float_ids = download_tech_files(float_ids)
%
% DESCRIPTION:
%   This function downloads the tech netcdf files for the floats with
%   the specified WMO IDs into subdirectory Tech. Extracting relevant
%   information from these files can be done outside the toolbox.
%
% INPUT:
%   float_ids : array with WMO IDs of the floats to be considered
%
% OUTPUT:
%   good_float_ids : WMO ID(s) of the float(s) whose tech files were downloaded
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

if nargin < 1
    warning('Usage: download_tech_files(float_ids)')
    return
end

% make sure Settings is initialized
if isempty(Settings)
    initialize_argo();
end

is_good = ones(length(float_ids), 1);
not_found = '';
count = 0;
for i = 1:length(float_ids)
    if ~download_float(float_ids(i), 'tech')
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
    fprintf('tech files could not be downloaded for floats:\n%s\n', ...
        not_found);
end
