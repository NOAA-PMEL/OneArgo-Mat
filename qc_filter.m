function [Data_good] = qc_filter(Data, variables, qc_flags, varargin)
% qc_filter  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   [Data_good] = qc_filter(Data, variables, qc_flags, varargin)
%
% DESCRIPTION:
%   This function generates a new data structure composed of chosen variables
%   based on provided QC flag values.
%
% INPUTS:
%   Data     : struct that must contain the given variables
%              (_ADJUSTED fields are used if available), as returned by
%              function load_float_data
%   variables: name(s) of the measured field(s)
%              This can be given as string (e.g., 'BBP700') for a single
%              variable or cell array for multiple variables, e.g.
%              {'DOXY';'NITRATE'}
%
% OPTIONAL INPUTS:
%   qc_flags: numerical array of QC flag values (default: [1,2])
%
%   Additional pairs of variables and qc flags can be specified, e.g.:
%     Data_good = qc_filter(Data, 'DOXY',[1:3,8], {'NITRATE';'CHLA'},[1,2]...
%         {'BBP532';'BBP470'},[1,2,8], 'DOWNWELLING_PAR',1);
%
% OUTPUT:
%   Data_good: struct that contains all the data from the input Data struct
%              that match the given QC flags
%
% AUTHORS:
%   N. Buzby (UW), H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL),
%   J. Plant, T. Maurer, Y. Takeshita (MBARI), D. Nicholson (WHOI),
%   and A. Gray (UW)
%
% CITATION:
%   H. Frenzel*, J. Sharp*, A. Fassbender, N. Buzby, J. Plant, T. Maurer,
%   Y. Takeshita, D. Nicholson, N. Buzby, A. Gray, 2021.
%   BGC-Argo-Mat: A MATLAB toolbox for accessing and
%   visualizing Biogeochemical Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.4971318.
%   (*These authors contributed equally to the code.)
%
% LICENSE: bgc_argo_mat_license.m
%
% DATE: December 1, 2021 (Version 1.1)

% assign default qc_flags if none provided as input
if nargin < 3
    qc_flags = [1,2];
end

if ischar(variables)
    variables = cellstr(variables);
end
nvar = length(variables);

% establish qc structure to reference
for v = 1:nvar
    qc_by_var.(variables{v}) = qc_flags;
end

% parse optional arguments
for i = 1:2:length(varargin)-1
    if ischar(varargin{i})
        qc_by_var.(varargin{i}) = varargin{i+1};
    else
        variables = varargin{i};
        nvar = length(variables);
        for v = 1:nvar
            qc_by_var.(variables{v}) = varargin{i+1};
        end
    end
end

variables = fields(qc_by_var);
nvar = length(variables);

floats = fieldnames(Data);
nfloats = length(floats);
for f = 1:nfloats
    % create basic structure to build off of
    Data_good.(floats{f}) = ...
        struct('CYCLE_NUMBER', Data.(floats{f}).CYCLE_NUMBER, ...
        'TIME', Data.(floats{f}).TIME, ...
        'LATITUDE', Data.(floats{f}).LATITUDE,...
        'LONGITUDE', Data.(floats{f}).LONGITUDE, ...
        'JULD', Data.(floats{f}).JULD);
    
    for v = 1:nvar
        if ~isfield(Data.(floats{f}),variables{v})
            warning('float %s does not contain variable %s', floats{f}, ...
                variables{v});
        else
            if isfield(Data.(floats{f}),[variables{v}, '_ADJUSTED']) && ...
                    sum(isfinite(Data.(floats{f}).([variables{v}, ...
                    '_ADJUSTED'])(:)))
                field = Data.(floats{f}).([variables{v}, '_ADJUSTED']);
                idx = ismember(Data.(floats{f}). ...
                    ([variables{v}, '_ADJUSTED_QC']),...
                    qc_by_var.(variables{v})) + 0; % convert to int
            else
                warning(['adjusted values for %s for are not available,',...
                    ' using raw values instead'], variables{v})
                field = Data.(floats{f}).(variables{v});
                idx = ismember(Data.(floats{f}).([variables{v}, '_QC']),...
                    qc_by_var.(variables{v})) + 0;
            end
            idx(~idx) = nan;
            Data_good.(floats{f}).(variables{v}) = field .* idx;
        end
    end
end
