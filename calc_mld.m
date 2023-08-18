function Data = calc_mld(Data, varargin)
% calc_mld  This function is part of the
% MATLAB toolbox for accessing Argo float Data.
%
% USAGE:
%   Data = calc_mld(Data, varargin)
%
% DESCRIPTION:
%   This function computes mixed layer depth using a density and/or
%   temperature threshold for all floats in the Data structure.
%
% INPUTS:
%   Data     : struct with data for one or more floats; expected fields are
%              PRES, TIME, LONGITUDE, LATITUDE, CYCLE_NUMBER, TEMP, PSAL
%              (or their _ADJUSTED equivalents)
%
% OPTIONAL INPUTS:
%   'raw',raw                     : use raw data if 'yes', adjusted data
%                                   if no (default: 'no')
%   'calc_mld_dens',calc_mld_dens : if set to 1, calculate mixed layer
%                                   depth (MLD) based on a density criterion
%   'dens_thres',dens_thres       : density threshold for MLD calculation;
%                                   default value is set in initialize_argo
%   'calc_mld_temp',calc_mld_temp : if set to 1, calculate mixed layer
%                                   depth based on a temperature criterion
%   'temp_thres',temp_thres       : temperature threshold for MLD calculation;
%                                   default value is set in initialize_argo
%   (Note that MLD can be computed both ways at the same time.)
%
% OUTPUT:
%   Data : struct with mixed layer depth added to the variables
%
% AUTHORS:
%   H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL), N. Buzby (UW)
%
% CITATION:
%   H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2022. OneArgo-Mat:
%   A MATLAB toolbox for accessing and visualizing Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.6588041
%
% LICENSE: oneargo_mat_license.m
%
% DATE: JUNE 1, 2022  (Version 1.0.1)

if nargin < 1
    warning('Usage: calc_mld(Data, [, varargin])')
    return
end

% set defaults
calc_mld_dens = 0;
calc_mld_temp = 0;
varargpass= {};
raw = 'no';

% parse optional arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'raw')
        raw = varargin{i+1};
    elseif strcmpi(varargin{i}, 'calc_mld_dens')
        calc_mld_dens = varargin{i+1};
    elseif strcmpi(varargin{i}, 'calc_mld_temp')
        calc_mld_temp = varargin{i+1};
    elseif strcmpi(varargin{i}, 'temp_thresh') || ...
            strcmpi(varargin{i}, 'dens_thresh')
        varargpass = [varargpass, varargin{i:i+1}];
    end
end

if ~calc_mld_dens && ~calc_mld_temp
    warning('Neither MLD method was selected!')
    return
end

floats = fieldnames(Data);
for f = 1:numel(floats)
    DataFloat = Data.(floats{f});
    fields = fieldnames(DataFloat);
    if calc_mld_dens
        if (strncmpi(raw,'y',1) && ~any(strcmp(fields, 'DENS'))) || ...
                (strncmpi(raw,'n',1) && ~any(strcmp(fields, 'DENS_ADJUSTED')))
            DataFloat = calc_auxil(DataFloat, 'calc_dens', 1, ...
                'raw', raw, 'float_name', floats{f});
        end
        DataFloat = calc_auxil(DataFloat, 'calc_mld_dens', 1, 'raw', raw, ...
            'float_name', floats{f}, varargpass);
    end
    if calc_mld_temp
        DataFloat = calc_auxil(DataFloat, 'calc_mld_temp', 1, 'raw', raw, ...
            'float_name', floats{f}, varargpass);
    end
    Data.(floats{f}) = DataFloat;
end
