function T = data2table(Data)
% data2table  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   T = data2table(Data)
%
% DESCRIPTION:
%   This function converts the Data structure to a flat table,
%   in which each row represents one observation of one profile.
%   The input Data struct can have one or more floats and any
%   number of variables. 
%   The output table will have columns for all variables that
%   occur at least in one of the floats. NaNs are used as fill values
%   for variables that some floats do not have.
%
% INPUT:
%   Data : struct (as returned from load_float_data)
%
% OUTPUT: 
%   T    : table (one observation per row, variables in columns)
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

% first step: determine all the fields across all floats
all_floats = fieldnames(Data);

all_fields = [];
for f = 1:length(all_floats)
    its_fields = fieldnames(Data.(all_floats{f}));
    all_fields = union(all_fields, its_fields);
end

T = table(); % this will be built up one float at a time

for f = 1:length(all_floats)
    [nlevels, nprof] = size(Data.(all_floats{f}).CYCLE_NUMBER);
    this_float = uint32(str2double(all_floats{f}(2:end)));
    Data.(all_floats{f}).WMOID = repmat(this_float, nlevels, nprof);
    its_fields = fieldnames(Data.(all_floats{f}));
    % these two fields have only one value per profile in the Data struct
    if any(strcmp(its_fields, 'MLD_DENS'))
        Data.(all_floats{f}).MLD_DENS = repmat(Data.(all_floats{f}).MLD_DENS,nlevels,1);
    end
    if any(strcmp(its_fields, 'MLD_TEMP'))
        Data.(all_floats{f}).MLD_TEMP = repmat(Data.(all_floats{f}).MLD_TEMP,nlevels,1);
    end
    for i = 1:length(its_fields)
        Data.(all_floats{f}).(its_fields{i}) = ...
            Data.(all_floats{f}).(its_fields{i})(:);
    end
    % all entries must have all fields; use NaNs as fill values
    non_fields = setdiff(all_fields,its_fields);
    for i = 1:length(non_fields)
        Data.(all_floats{f}).(non_fields{i}) = ...
            nan(size( Data.(all_floats{f}).TIME));
    end
    T = cat(1,T,struct2table(Data.(all_floats{f})));
end

