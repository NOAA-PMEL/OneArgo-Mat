function all_vars = combine_variables(base_vars, new_vars)
% combine_variables  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   all_vars = combine_variables(base_vars, new_vars)
%
% DESCRIPTION:
%   This function combines the given variables along with all
%   associated variables (e.g., _ADJUSTED, _QC, etc.) and returns them.
%
% INPUTS:
%   base_vars : cell array with names of the basic variables
%   new_vars  : cell array with names of additional variables (tracer
%               fields etc.) that have the standard associated variables
%
% OUTPUT:
%  all_vars   : cell array with names of all variables
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

cnt_vars = length(base_vars);

nvars = cnt_vars + 6*(length(new_vars)) - sum(ismember('PRES', new_vars));
all_vars = cell(nvars, 1); % pre-allocated
all_vars(1:cnt_vars) = base_vars;

% include associated variables
for i = 1:length(new_vars)
    cnt_vars = cnt_vars+1;
    all_vars{cnt_vars} = new_vars{i};
    cnt_vars = cnt_vars+1;
    all_vars{cnt_vars} = [new_vars{i}, '_QC'];
    if ~strcmp(new_vars{i}, 'PRES')
        cnt_vars = cnt_vars+1;
        all_vars{cnt_vars} = [new_vars{i}, '_dPRES'];
    end
    cnt_vars = cnt_vars+1;
    all_vars{cnt_vars} = [new_vars{i}, '_ADJUSTED'];
    cnt_vars = cnt_vars+1;
    all_vars{cnt_vars} = [new_vars{i}, '_ADJUSTED_QC'];
    cnt_vars = cnt_vars+1;
    all_vars{cnt_vars} = [new_vars{i}, '_ADJUSTED_ERROR'];
    cnt_vars = cnt_vars+1;
    all_vars{cnt_vars} = ['PROFILE_', new_vars{i}, '_QC'];
end
