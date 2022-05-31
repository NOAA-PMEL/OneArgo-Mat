function [n_prof, n_param, n_levels] = get_dims(filename)
% get_dims  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   [n_prof, n_param, n_levels] = get_dims(filename)
%
% DESCRIPTION:
%   This function determines the number of profiles, parameters,
%   and depth levels in an Sprof netcdf file.
%
% INPUT:
%   filename    : the name of the Sprof file
%
% OUTPUTS:
%   n_prof      : number of profiles (integer)
%   n_param     : number of parameters (integer)
%   n_levels    : number of depth levels (integer)
%
% AUTHORS:
%   J. Sharp, H. Frenzel, A. Fassbender (NOAA-PMEL), N. Buzby (UW)
%
% CITATION:
%   H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2022. OneArgo-Mat:
%   A MATLAB toolbox for accessing and visualizing Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.6588042
%
% LICENSE: oneargo_mat_license.m
%
% DATE: JUNE 1, 2022  (Version 1.0.1)

info = ncinfo(filename); % Read netcdf information
dims = info.Dimensions; % Extract dimensional information

% Determine names of dimensional properties
dimensions = {dims.('Name')};

% Find 'number of profiles', 'number of parameters',
% and 'number of depth levels'
n_prof = dims(strcmp(dimensions, 'N_PROF')).Length;
n_param = dims(strcmp(dimensions, 'N_PARAM')).Length;
n_levels = dims(strcmp(dimensions, 'N_LEVELS')).Length;



