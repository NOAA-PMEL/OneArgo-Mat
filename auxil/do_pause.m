function do_pause()
% do_pause  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   do_pause()
%
% DESCRIPTION:
%   It asks the user to hit ENTER if Settings.use_pause is non-zero.
%   Otherwise, flow control returns the caller.
%
% AUTHORS: 
%   H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL),
%   J. Plant, T. Maurer, Y. Takeshita (MBARI), D. Nicholson (WHOI),
%   and A. Gray (UW)
%
% CITATION:
%   BGC-Argo-Mat: A MATLAB toolbox for accessing and visualizing
%   Biogeochemical Argo data,
%   H. Frenzel*, J. Sharp*, A. Fassbender, J. Plant, T. Maurer, 
%   Y. Takeshita, D. Nicholson, and A. Gray; 2021
%   (*These authors contributed equally to the code.)
%
% LICENSE: bgc_argo_mat_license.m
%
% DATE: June 15, 2021

global Settings;

if Settings.use_pause
    disp('Please hit ENTER to continue');
    pause;
end
