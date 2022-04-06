function initialize_tech(file_name)
% initialize_tech  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
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
%   H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL), N. Buzby (UW),
%   J. Plant, T. Maurer, Y. Takeshita (MBARI), D. Nicholson (WHOI),
%   and A. Gray (UW)
%
% CITATION:
%   H. Frenzel*, J. Sharp*, A. Fassbender, N. Buzby, J. Plant, T. Maurer,
%   Y. Takeshita, D. Nicholson, A. Gray, 2021. BGC-Argo-Mat: A MATLAB
%   toolbox for accessing and visualizing Biogeochemical Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.4971318.
%   (*These authors contributed equally to the code.)
%
% LICENSE: bgc_argo_mat_license.m
%
% DATE: FEBRUARY 22, 2022  (Version 1.2)

global Tech;

fid = fopen(file_name);
H = textscan(fid,'%s %s %s','headerlines',9,...
    'delimiter',',','whitespace','');
fclose(fid);
Tech.file_path = H{1};
tech_wmoid = regexp(Tech.file_path,'\d{7}','once','match');
Tech.file_name = regexprep(tech_wmoid,'\d{7}','$0_tech.nc');
Tech.update = H{3};
Tech.wmoid = str2double(tech_wmoid);
