function initialize_meta(file_name)
% initialize_meta  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   initialize_meta(file_name)
%
% DESCRIPTION:
%   This function initializes the global struct Meta by reading
%   the index file and processing its information.
%
% INPUTS:
%   file_name : name of the index file (with local path)
%
% OUTPUTS: None. Meta is filled with fields.
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

global Meta;

fid = fopen(file_name);
H = textscan(fid,'%s %s %s %s','headerlines',9,...
    'delimiter',',','whitespace','');
fclose(fid);
Meta.file_path = H{1};
meta_wmoid = regexp(Meta.file_path,'\d+','once','match');
Meta.file_name = regexprep(Meta.file_path, '[a-z]+/\d+/', '');
Meta.update = H{4};
Meta.wmoid = str2double(meta_wmoid);
