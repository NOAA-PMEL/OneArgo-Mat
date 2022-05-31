function initialize_traj(file_name)
% initialize_traj  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   initialize_traj(file_name)
%
% DESCRIPTION:
%   This function initializes the global struct Traj by reading
%   the index file and processing its information.
%
% INPUTS:
%   file_name : name of the index file (with local path)
%
% OUTPUTS: None. Traj is filled with fields.
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

global Traj;

fid = fopen(file_name);
H = textscan(fid,'%s %f %f %f %f %s %s %s','headerlines',9,...
    'delimiter',',','whitespace','');
fclose(fid);
Traj.file_path = H{1};
Traj.lat_max = H{2};
Traj.lat_min = H{3};
Traj.lon_max = H{4};
Traj.lon_min = H{5};
Traj.profiler_type = H{6};
traj_wmoid = regexp(Traj.file_path,'\d+','once','match');
Traj.file_name = regexprep(Traj.file_path, '[a-z]+/\d+/', '');
Traj.update = H{8};
Traj.wmoid = str2double(traj_wmoid);
