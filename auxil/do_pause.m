function do_pause()
% do_pause  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   do_pause()
%
% DESCRIPTION:
%   It asks the user to hit ENTER if Settings.use_pause is non-zero.
%   Otherwise, flow control returns to the caller.
%
% INPUT: None.
%
% OUTPUT: None.
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

% make sure Settings is initialized
if isempty(Settings)
    initialize_argo();
end

if Settings.use_pause
    disp('Please hit ENTER to continue');
    pause;
end
