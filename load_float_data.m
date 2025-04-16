function [Data, Mdata] = load_float_data(float_ids, variables, ...
    float_profs, varargin)
% load_float_data  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   [Data, Mdata] = load_float_data(float_ids [, variables] ...
%                                   [, float_profs], varargin)
%
% DESCRIPTION:
%   This function loads data (at least one variable)
%   of at least one specified float.
%
% INPUT:
%   float_ids   : WMO ID(s) of one or more floats
%
% OPTIONAL INPUTS:
%   variables   : cell array with variable names to be loaded (use 'ALL'
%                 to load all available variables, which may differ by
%                 float)
%   float_profs : cell array with indices of selected profiles (per float,
%                 not global)
%   'interp_lonlat', intp : if intp is 'yes' (default), missing lon/lat
%                 values (e.g., under ice) will be interpolated;
%                 set intp to 'no' to suppress interpolation;
%                 the default is taken from Settings.interp_lonlat (defined
%                 in initialize_argo.m)
%
% OUTPUTS:
%   Data        : struct with the requested variables (including QC flags.
%                 adjusted values if available) and general ones
%                 (LONGITUDE, LATITUDE, JULD, etc.)
%   Mdata       : struct with meta data (WMO_NUMBER, PI_NAME, etc.)
%
% AUTHORS:
%   J. Sharp and H. Frenzel (UW-CICOES), A. Fassbender (NOAA-PMEL), N. Buzby (UW)
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

add_pres = 0; % default: do not add 'PRES' to list of variables

if nargin < 1
    warning('Usage: load_float_data(float_ids [, variables, float_profs])')
    return
end
if nargin >= 2
    % convert requested variable to cell array if necessary (string was used)
    if ischar(variables)
        variables = cellstr(variables);
    end
    if ~sum(ismember('PRES', variables))
        add_pres = 1;
    end
else
    variables = {};
end
if nargin < 3
    % by default, all profiles of the given floats are loaded
    float_profs = [];
end

interp_ll = Settings.interp_lonlat;
% parse optional arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'interp_lonlat')
        interp_ll = varargin{i+1};
    else
        warning('unknown option: %s', varargin{i+1});
    end
end

% INITIALIZE STRUCTURES FOR OUTPUT
Data = struct();
Mdata = struct();
fields_mdata = {'PROJECT_NAME';'PI_NAME';'DATA_CENTRE'};

% only some variables are always loaded, others only by request
base_vars = {'CYCLE_NUMBER'; 'DIRECTION'; 'JULD'; 'JULD_QC'; ...
    'JULD_LOCATION'; 'LATITUDE'; 'LONGITUDE'; 'POSITION_QC'; ...
    'PARAMETER_DATA_MODE'; 'PARAMETER'; 'DATA_MODE'};

if ~isempty(variables) && strcmp(variables{1}, 'ALL')
    use_all_vars = 1;
    base_vars{end+1} = 'PROFILE_PRES_QC';
else
    use_all_vars = 0;
    if add_pres
        % if no variables are specified (e.g., to plot trajectories only),
        % loading pressure and associated variables is not needed
        variables{end+1} = 'PRES';
        base_vars{end+1} = 'PROFILE_PRES_QC';
    end
    
    add_vars = ismember(Settings.avail_vars, variables);
    new_vars = Settings.avail_vars(add_vars);
    all_vars = combine_variables(base_vars, new_vars);
end

% download Sprof files if necessary
good_float_ids = download_multi_floats(float_ids);

if length(good_float_ids) < length(float_ids) && ~isempty(float_profs)
    [~, idx] = intersect(float_ids, good_float_ids, 'stable');
    float_profs = float_profs(idx);
end

% LOOP TO IMPORT PROFILES AND EXTRACT VARIABLES
for n = 1:length(good_float_ids)
    floatnum = good_float_ids(n);
    is_bgc_float = strcmp(Float.type{Float.wmoid == floatnum}, 'bgc');
    str_floatnum = ['F', num2str(floatnum)];
    filename = sprintf('%s%s', Settings.prof_dir, ...
        Float.file_name{Float.wmoid == floatnum});
    if use_all_vars
        info = ncinfo(filename); % Read netcdf information
        these_vars = {info.Variables.('Name')};
        add_vars = ismember(Settings.avail_vars, these_vars);
        new_vars = Settings.avail_vars(add_vars);
        all_vars = combine_variables(base_vars, new_vars);
    end
    % Find 'number of profiles', 'number of parameters', and 'number of
    % depth levels'
    [n_prof, n_param, n_levels] = get_dims(filename);
    n_vars = length(all_vars);
    % Extract data from netcdf file and save data in proper structures
    for l = 1:n_vars
        if ~is_bgc_float && (strcmp(all_vars{l}, 'PARAMETER_DATA_MODE') ...
                || endsWith(all_vars{l}, '_dPRES'))
            continue; % these variables are only in Sprof files
        end
        if is_bgc_float && strcmp(all_vars{l}, 'DATA_MODE')
            continue; % this variable is only in prof files
        end
        try
            tmp = ncread(filename,all_vars{l});
        catch
            warning('%s not found in %s', all_vars{l}, filename)
            continue; % skip this variable for this float
        end
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
        if isequal(size(tmp), [n_levels, n_prof]) && ...
                ~strcmp(all_vars{l}, 'DATA_MODE')
            Data.(str_floatnum).(all_vars{l}) = tmp;
        elseif isequal(size(tmp), [n_prof 1])
            if strcmp(all_vars{l}, 'DATA_MODE')
                % emulate the setup for bgc floats
                Mdata.(str_floatnum).PARAMETER_DATA_MODE = ...
                    repmat(tmp',[n_param, 1]);
            else
                Data.(str_floatnum).(all_vars{l}) = repmat(tmp', n_levels, 1);
            end
        else
            chars = permute(sum(sum(tmp)), [3 4 1 2]);
            [~, idx] = max(chars(end,:));
            Mdata.(str_floatnum).(all_vars{l}) = tmp(:,:,1,idx);
        end
        clear tmp;
    end
    if strncmpi(interp_ll, 'y', 1)
        Data = interp_lonlat(Data, floatnum);
    end
    
    % Add WMO float number to metadata
    Mdata.(str_floatnum).WMO_NUMBER = floatnum;
    
    % parse parameter names
    temp = cell(n_param, 1);
    % extract parameter names as coherent strings
    for m = 1:n_param
        temp{m} = ...
            strrep(Mdata.(str_floatnum).('PARAMETER')(:,m)', ' ', '');
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
    Mdata.(str_floatnum) = rmfield(Mdata.(str_floatnum),...
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
