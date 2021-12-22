function good_variables = check_variables(variables)
% check_variables  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   good_variables = check_variables(variables)
%
% DESCRIPTION:
%   This function checks if the specified variable(s) are available.
%   The available variables are returned.
%
% INPUTS:
%   variables : cell array of variable(s) (or one variable as a string)
%
% OUTPUTS:
%   good_variables : cell array of available variable(s)
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
% DATE: DECEMBER 1, 2021  (Version 1.1)

global Settings;

% make sure Settings is initialized
if isempty(Settings)
    initialize_argo();
end

if ischar(variables)
    variables = cellstr(variables);
end

good_variables = variables;

for i = 1:length(variables)
    if ~any(strcmp(variables{i}, Settings.avail_vars))
        good_variables(strcmp(variables{i}, good_variables)) = [];
    end
end