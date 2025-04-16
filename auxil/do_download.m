function tf = do_download(dest_path)
% do_download  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   tf = do_download(dest_path)
%
% DESCRIPTION:
%   This function determines if a file should be downloaded or not
%   (i.e., if it exists already at the given dest_path), based on the
%   "update" option from the global Settings structure.
%
% INPUT:
%   dest_path : local destination path for a file, which may not yet exist
%
% OUTPUT:
%   tf        : True (1) or false (0) - is download needed? (This depends
%               on the value of Settings.update and the existence and age
%               of the file.)
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

global Settings;

if nargin < 1
    warning('Usage: do_download(dest_path)')
    tf = 1;
    return
end

if exist(dest_path, 'file') ~= 2
    tf = 1;
elseif Settings.update == 0 || Settings.update == 1
    tf = Settings.update;
else
    file_info = dir(dest_path);
    file_age = (now - datenum(file_info.date))*86400;
    tf = file_age > Settings.update;
end
