function dacs = list_dacs(float_ids, mode)
% list_dacs  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   dacs = list_dacs(float_ids [, mode])
%
% DESCRIPTION:
%   This function lists the DACs (Data Assembly Centers) that handle
%   the data of the given floats.
%
% INPUT:
%   float_ids : array with WMO IDs of the floats to be considered
%
% OPTIONAL INPUT:
%   mode : either 'summary' (default; return a cell array with
%          all DACs handling the specified floats)
%          or 'all' (also displays a table with a listing of
%          the DAC for all specified floats)
%
% OUTPUT:
%   dacs : cell array with DACs
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

global Float Settings;

% make sure Settings is initialized
if isempty(Settings)
    initialize_argo();
end

if nargin < 2
    mode = 'summary';
end

fidx = arrayfun(@(x) find(Float.wmoid==x, 1), float_ids, ...
    'UniformOutput', false);
fidx(cellfun(@isempty, fidx)) = [];
fidx = cell2mat(fidx);
if isempty(fidx)
    warning('no valid float IDs specified')
    dacs = {};
    return
end

all_dacs = Float.dac(fidx);
dacs = unique(all_dacs);

if strcmpi(mode, 'all')
    fprintf('WMOID     DAC\n')
    fprintf('%s\n', repelem('-', 20));
    for f = 1:length(fidx)
        fprintf('%-7d   %s\n', Float.wmoid(fidx(f)), Float.dac{fidx(f)});
    end
end