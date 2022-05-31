function num_good_values = plot_one_profile(axes, data, pres, data_qc, ...
    qc_flags, plot_obs, plot_color)
% plot_one_profile  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   num_good_values = plot_one_profile(axes, data, pres, data_qc, ...
%       qc_flags, plot_obs [, plot_color])
%
% DESCRIPTION:
%   This function plots one profile of one float for one variable.
%
% INPUTS:
%   axes     : axes of the plot window, which must exist already
%   data     : values of the variable (vector)
%   pres     : pressure values (vector, same size as data)
%   data_qc  : QC values of the data (vector, same size as data)
%   qc_flags : the allowed QC flags (only data whose data_qc values
%              match the specified qc_flags will be shown)
%   plot_obs : switch to plot markers at depths of observations (1);
%              or not (0)
%
% OPTIONAL INPUT:
%   plot_color : color to use for the line and markers (default: black);
%              standard Matlab short format colors can be used (e.g.: 'r')
%
% OUTPUT:
%   num_good_values : number of values (data and pres) that are finite
%              and whose data_qc values match specified qc_flags
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

if nargin < 7
    plot_color = 'k';
end

idx = isfinite(data) & isfinite(pres) & ismember(data_qc, qc_flags);
if sum(idx)
    plot(axes, data(idx), pres(idx), 'color', plot_color);
    if plot_obs
        scatter(axes, data(idx), pres(idx), 2, plot_color, 'filled');
    end
end
num_good_values = sum(idx);
