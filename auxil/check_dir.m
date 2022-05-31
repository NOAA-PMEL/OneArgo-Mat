function success = check_dir(ddir)
% check_dir  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   success = check_dir(ddir)
%
% DESCRIPTION:
%   This function determines if a directory needs to be created and
%   does so if necessary.
%
% INPUT:
%   ddir : directory (this can be a relative or absolute path)
%
% OUTPUT:
%   success : 0 if ddir did not exist yet and cannot be created; 1
%             otherwise
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
    warning('Usage: check_dir(ddir)')
    return
end

if ~exist(ddir, 'dir') && ~mkdir(ddir)
    warning('Could not create directory %s', ddir)
    success = 0;
else
    success = 1;
end
