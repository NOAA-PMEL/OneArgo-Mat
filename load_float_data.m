function [Data, Mdata] = load_float_data(float_ids, variables, float_profs)
% load_floats  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   [Data, Mdata] = load_float_data(float_ids, variables, float_profs)
%
% DESCRIPTION:
%   This function loads data (at least one variable)
%   of at least one specified float.
%
% INPUTS:
%   float_ids   : WMO ID of one or more floats
%
% OPTIONAL INPUTS:
%   variables   : cell array with variable names to be loaded
%   float_profs : cell array with IDs of selected profiles (per float,
%                 not global)
%
% OUTPUT:
%   Data        : struct with the requested variables (including QC flags.
%                 adjusted values if available) and general ones
%                 (LONGITUDE,LATITUDE,JULD etc.)
%   Mdata       : struct with meta data (WMO_NUMBER)
%
% AUTHORS: 
%   J. Sharp, H. Frenzel, A. Fassbender (NOAA-PMEL),
%   J. Plant, T. Maurer, Y. Takeshita (MBARI), D. Nicholson (WHOI),
%   and A. Gray (UW)
%
% CITATION:
%   H. Frenzel*, J. Sharp*, A. Fassbender, J. Plant, T. Maurer,
%   Y. Takeshita, D. Nicholson, A. Gray, 2021. BGC-Argo-Mat: A MATLAB
%   toolbox for accessing and visualizing Biogeochemical Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.4971318.
%   (*These authors contributed equally to the code.)
%
% LICENSE: bgc_argo_mat_license.m
%
% DATE: October 15, 2021 (Version 1.1)

global Settings;

if nargin < 1
    warning('Usage: load_float_data(float_ids [, variables, float_profs])')
    return
end
if nargin >= 2
    % convert requested variable to cell array if necessary (string was used)
    if ischar(variables)
        variables = cellstr(variables);
    end
    add_pres = 1;
else
    variables = {};
    add_pres = 0;
end
if nargin < 3
    % by default, all profiles of the given floats are loaded
    float_profs = [];
end

% INITIALIZE STRUCTURES FOR OUTPUT
Data = struct();
Mdata = struct();

% pre-allocate all_vars cell array
nvars = 9 + 6*(length(variables) + add_pres);
all_vars = cell(nvars, 1);

% only some variables are always loaded, others only by request
all_vars(1:9) = {'CYCLE_NUMBER'; 'DIRECTION'; 'JULD'; 'JULD_QC'; ...
    'JULD_LOCATION'; 'LATITUDE'; 'LONGITUDE'; 'PARAMETER_DATA_MODE'; ...
    'PARAMETER'};

if add_pres
    % if no variables are specified (e.g., to plot trajectories only),
    % loading pressure and associated variables is not needed
    variables{end+1} = 'PRES';
    all_vars{10} = 'PROFILE_PRES_QC';
end
cnt_vars = 9 + add_pres;

fields_mdata = {'PROJECT_NAME';'PI_NAME';'DATA_CENTRE'};

add_vars = ismember(Settings.avail_vars, variables);
new_vars = Settings.avail_vars(add_vars);

% always include all associated variables
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
end

% download Sprof files if necessary
good_float_ids = download_multi_floats(float_ids);

% LOOP TO IMPORT PROFILES AND EXTRACT VARIABLES
for n = 1:length(good_float_ids)
    floatnum = good_float_ids(n);
    str_floatnum = ['F', num2str(floatnum)];
    filename = sprintf('%s%d_Sprof.nc', Settings.prof_dir, floatnum);
    
    % Find 'number of profiles', 'number of parameters', and 'number of
    % depth levels'
    [n_prof, n_param, n_levels] = get_dims(filename);
    n_vars = length(all_vars);
    % Extract data from netcdf file and save data in proper structures
    for l = 1:n_vars
        tmp = ncread(filename,all_vars{l});
        % CONVERT QUALITY FLAGS TO NUMERIC FORMAT
        if endsWith(all_vars{l},'_QC') && ... % Check for QC identifier
                ~startsWith(all_vars{l},'PROFILE') % But not a profile QC
            if isequal(size(tmp), [n_prof 1])
                tmp = repmat(tmp',n_levels,1);
            end
            tmp = tmp(:);
            tmp = strrep(tmp', ' ', '0')';
            tmp = str2num(tmp);
            tmp = reshape(tmp, n_levels, n_prof);
        end
        if isequal(size(tmp), [n_levels, n_prof])
            Data.(str_floatnum).(all_vars{l}) = tmp;
        elseif isequal(size(tmp), [n_prof 1])
            Data.(str_floatnum).(all_vars{l}) = repmat(tmp', n_levels, 1);
        else
            Mdata.(str_floatnum).(all_vars{l}) = tmp;          
        end
        clear tmp;
    end
        
    % Add WMO float number to metadata
    Mdata.(str_floatnum).WMO_NUMBER = floatnum;

    % parse parameter names
    temp = cell(n_param, 1);
    % extract parameter names as coherent strings
    for m = 1:n_param
        temp{m} = ...
            strrep(Mdata.(str_floatnum).('PARAMETER')(:,m,1,1)',' ','');
    end
    params_keep = ismember(temp,new_vars);
    Mdata.(str_floatnum).('PARAMETER') = temp(params_keep);
    clear temp;
    
    % parse parameter data modes
    % create data mode variable for each parameter
    % expand that variable to match size of data matrix
    p = 0; % index of parameters actually used
    for m = 1:n_param
        if params_keep(m)
            p = p + 1;
            Data.(str_floatnum).(...
                [cell2mat(Mdata.(str_floatnum).PARAMETER(p)),...
                '_DATA_MODE']) = ...
                repmat(Mdata.(str_floatnum).('PARAMETER_DATA_MODE')(m,:),...
                n_levels,1);
        end
    end
    
    % clear both parameter and parameter data mode from metadata
    Mdata.(str_floatnum) = ...
        rmfield(Mdata.(str_floatnum),...
        {'PARAMETER','PARAMETER_DATA_MODE'});
    
    % add information about deploying organization and PI to meta data
    for f = 1:length(fields_mdata)
        this_field = ncread(filename, fields_mdata{f});
        Mdata.(str_floatnum).(fields_mdata{f}) = ...
            strcat(this_field(:,end)');
    end

    % CONVERT JULD VARIABLE TO SERIAL DATE (SINCE YEAR 1950)
    % AND SAVE AS 'TIME'
    Data.(str_floatnum).('TIME') = ...
        datenum(Data.(str_floatnum).('JULD'))+...
        datenum([1950 1 1]);
    
    % Select only specified profiles
    if ~isempty(float_profs)
        names = fieldnames(Data.(str_floatnum));
        for l = 1:numel(names)
            Data.(str_floatnum).(names{l}) = ...
                Data.(str_floatnum).(names{l})(:,float_profs{n});
        end
    end
end

