function initialize_tech(file_name)
% initialize_tech  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   initialize_tech(file_name)
%
% DESCRIPTION:
%   This function initializes the global struct Tech by reading
%   the index file and processing its information.
%
% INPUTS:
%   file_name : name of the index file (with local path)
%
% OUTPUTS: None. Tech is filled with fields.
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

global Tech;

fid = fopen(file_name);
H = textscan(fid,'%s %s %s','headerlines',9,...
    'delimiter',',','whitespace','');
fclose(fid);
Tech.file_path = H{1};
tech_wmoid = regexp(Tech.file_path,'\d+','once','match');
Tech.file_name = regexprep(Tech.file_path, '[a-z]+/\d+/', '');
Tech.update = H{3};
Tech.wmoid = str2double(tech_wmoid);
