function Data_good = qc_filter(Data, variables, qc_flags, varargin)
% qc_filter  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   Data_good = qc_filter(Data, variables [, qc_flags], varargin)
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
%   'raw',raw : if set to 'n' (default), use only adjusted data; however,
%             if not all floats have adjusted data for the selected
%             variables, raw data will be used;
%             if set to 'y', only raw data will be used;
%             if set to 'no_strict', only adjusted data will be used;
%             floats without adjusted data for the selected variables
%             will be omitted from the returned struct
%
% OUTPUT:
%   Data_good: struct that contains all the variables from the input Data
%              struct that match the given QC flags;
%              all other values are set to NaN (the size of the arrays is
%              unchanged)
%
% AUTHORS:
%   N. Buzby (UW), H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL)
%
% CITATION:
%   H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2022. OneArgo-Mat:
%   A MATLAB toolbox for accessing and visualizing Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.6588041
%
% LICENSE: oneargo_mat_license.m
%
% DATE: JUNE 1, 2022  (Version 1.0.1)

% assign default qc_flags if none provided as input
if nargin < 3
    qc_flags = [1,2];
end
raw = 'n'; % default: do not take raw values if adjusted are available

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
    if strcmpi(varargin{i}, 'raw')
        raw = varargin{i+1};
    elseif ischar(varargin{i})
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
        'POSITION_QC', Data.(floats{f}).POSITION_QC, ...
        'JULD', Data.(floats{f}).JULD);
    
    for v = 1:nvar
        if ~isfield(Data.(floats{f}),variables{v})
            warning('float %s does not contain variable %s', floats{f}, ...
                variables{v});
        else
            if strncmpi(raw, 'y', 1)
                field = Data.(floats{f}).(variables{v});
                idx = ismember(Data.(floats{f}).([variables{v}, '_QC']),...
                    qc_by_var.(variables{v})) + 0;
            elseif isfield(Data.(floats{f}),[variables{v}, '_ADJUSTED']) && ...
                    sum(isfinite(Data.(floats{f}).([variables{v}, ...
                    '_ADJUSTED'])(:)))
                field = Data.(floats{f}).([variables{v}, '_ADJUSTED']);
                idx = ismember(Data.(floats{f}). ...
                    ([variables{v}, '_ADJUSTED_QC']),...
                    qc_by_var.(variables{v})) + 0; % convert to int
            elseif strcmpi(raw, 'no_strict')
                % if one variable doesn't have any good values, do not
                % use any data for this float
                warning('for float %s, adjusted values for %s are not available,%s',...
                    floats{f}, variables{v}, ' this float will not be used');
                Data_good = rmfield(Data_good, floats{f});
                break;
            elseif strncmpi(raw, 'n', 1)
                warning('for float %s, adjusted values for %s are not available,%s',...
                    floats{f}, variables{v}, ' this field will not be used');
                continue;
            else
                warning(['adjusted values for %s for are not available,',...
                    ' using raw values instead'], variables{v})
                field = Data.(floats{f}).(variables{v});
                idx = ismember(Data.(floats{f}).([variables{v}, '_QC']),...
                    qc_by_var.(variables{v})) + 0;
            end
            idx(~idx) = nan;
            if strncmpi(raw, 'y', 1)
                Data_good.(floats{f}).(variables{v}) = field .* idx;
                Data_good.(floats{f}).([variables{v}, '_QC']) = ...
                    Data.(floats{f}).([variables{v}, '_QC']);
            else % adjusted values were requested
                Data_good.(floats{f}).([variables{v}, '_ADJUSTED']) = ...
                    field .* idx;
                Data_good.(floats{f}).([variables{v}, '_ADJUSTED_QC']) = ...
                    Data.(floats{f}).([variables{v}, '_ADJUSTED_QC']);
            end
        end
    end
end
