function list_sensors()
% list_sensors  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   list_sensors()
%
% DESCRIPTION:
%   This function lists all currently defined sensors by their
%   short names (e.g., DOXY) and long names (e.g., Oxygen) along with
%   the units.
%
% INPUTS:
%   None.
%
% OUTPUTS: 
%   No return values. Short and long names as well as units are displayed
%   in the command window.
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

nvars = length(Settings.avail_vars);

fprintf('\nThese are the currently defined sensors:\n\n')
fprintf('%-20s   %-42s  %s\n', 'Short name', 'Long name', 'Units');
fprintf('%s\n', repelem('-', 80));
for i = 1:nvars
    [long_name, units] = get_var_name_units(Settings.avail_vars{i});
    fprintf('%-20s   %-42s  %s\n', Settings.avail_vars{i}, long_name, ...
        strrep(strrep(strrep(units, '{', ''), '}', ''), '\mu', 'u'))
end
