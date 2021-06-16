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
%                 (if not set: a default float is used as a demo)
%
% OPTIONAL INPUTS:
%   variables   : cell array with variable names to be loaded
%   float_profs : cell array with IDs of selected profiles (per float,
%                 not global)
%
% OUTPUT:
%   Data        : struct with the requested variables (including QC flags.
%                 adjusted values if available) and general ones
%                 (LONGITUDE,LATITUDE,JULD)
%   Mdata       : struct with meta data (WMO_NUMBER)
%
% AUTHORS: 
%   J. Sharp, H. Frenzel, A. Fassbender (NOAA-PMEL),
%   J. Plant, T. Maurer, Y. Takeshita (MBARI), D. Nicholson (WHOI),
%   and A. Gray (UW)
%
% CITATION:
%   BGC-Argo-Mat: A MATLAB toolbox for accessing and visualizing
%   Biogeochemical Argo data,
%   H. Frenzel*, J. Sharp*, A. Fassbender, J. Plant, T. Maurer, 
%   Y. Takeshita, D. Nicholson, and A. Gray; 2021
%   (*These authors contributed equally to the code.)
%
% LICENSE: bgc_argo_mat_license.m
%
% DATE: June 15, 2021

global Settings;

if nargin < 1
    warning('Usage: load_float_data(float_ids [, variables, float_profs])')
end

% only some variables are always loaded, others only by request
all_vars = {'CYCLE_NUMBER'; 'DIRECTION'; 'JULD'; 'JULD_QC'; ...
    'JULD_LOCATION'; 'LATITUDE'; 'LONGITUDE'; 'PARAMETER_DATA_MODE'; ...
    'PARAMETER'};

if nargin >= 2
    % convert requested variable to cell array if necessary (string was used)
    if ischar(variables)
        variables = cellstr(variables);
    end
    % if no variables are specified (e.g., to plot trajectories),
    % loading pressure and associated variables is not needed
    variables{end+1} = 'PRES';
    all_vars{end+1} = 'PROFILE_PRES_QC';
else
    variables = {};
end
if nargin < 3
    % by default, all profiles of the given floats are loaded
    float_profs = [];
end

% INITIALIZE STRUCTURES FOR Data OUTPUT
Data = struct();
Mdata = struct();

add_vars = ismember(Settings.avail_vars, variables);
new_vars = Settings.avail_vars(add_vars);

% always include all associated variables
for i = 1:length(new_vars)
    all_vars{end+1} = new_vars{i};
    all_vars{end+1} = [new_vars{i}, '_QC'];
    if ~strcmp(new_vars{i}, 'PRES')
        all_vars{end+1} = [new_vars{i}, '_dPRES'];
    end
    all_vars{end+1} = [new_vars{i}, '_ADJUSTED'];
    all_vars{end+1} = [new_vars{i}, '_ADJUSTED_QC'];
    all_vars{end+1} = [new_vars{i}, '_ADJUSTED_ERROR'];
end

% download Sprof files if necessary
good_float_ids = download_multi_floats(float_ids);

% LOOP TO IMPORT PROFILES AND EXTRACT VARIABLES
for n = 1:length(good_float_ids)
    floatnum = good_float_ids(n);
    filename = sprintf('%s%d_Sprof.nc', Settings.prof_dir, floatnum);
    
    % LOAD VARIABLES FROM FILE
    info = ncinfo(filename); % Read netcdf information
    dims = info.Dimensions; % Extract dimensional information
    % Determine names of dimensional properties
    dimensions = cell(numel(dims),1);
    for h=1:numel(dims)
        dimensions(h) = {dims(h).Name};
    end
    % Find 'number of profiles', 'number of parameters', and 'number of
    % depth levels'
    profidx  = find(strcmp(dimensions,'N_PROF'));
    n_prof = dims(profidx).Length;
    paramidx = find(strcmp(dimensions,'N_PARAM'));
    n_param = dims(paramidx).Length;
    levidx   = find(strcmp(dimensions,'N_LEVELS'));
    n_levels = dims(levidx).Length;
    amt = length(all_vars);
    names = cell(amt,1); % Pre-allocate variable names
    mnames = cell(amt,1); % Pre-allocate meta-variable names
    for l=1:numel(names)
        names(l) = cellstr(all_vars{l});
        mnames(l) = names(l);
    end
    % Extract data from netcdf, log variable names, and save data in
    % proper structures
    for l=1:numel(names)
        % Read in data
        Data.(strcat('F',num2str(floatnum))).(names{l}) = ...
            ncread(filename,char(names(l)));
        Mdata.(strcat('F',num2str(floatnum))).(mnames{l}) = ...
            Data.(strcat('F',num2str(floatnum))).(names{l});
        % For measured variables
        if numel(size(Data.(strcat('F',...
                num2str(floatnum))).(names{l}))) == 2 && ...
                all(size(Data.(strcat('F',...
                num2str(floatnum))).(names{l})) == [n_levels n_prof])
            % Remove metadata fields
            Mdata.(strcat('F',num2str(floatnum))) = ...
                rmfield(Mdata.(strcat('F',num2str(floatnum))),mnames{l});
            mnames{l} = [];
            % For descriptive meta variables (1 value per profile)
        elseif numel(size(Data.(strcat('F',...
                num2str(floatnum))).(names{l}))) == 2 && ...
                all(size(Data.(strcat('F',...
                num2str(floatnum))).(names{l})) == [n_prof 1])
            % Replicate to match number of observations
            Data.(strcat('F',num2str(floatnum))).(names{l}) = ...
                repmat(Data.(strcat('F',...
                num2str(floatnum))).(names{l})',n_levels,1);
            % Remove metadata fields
            Mdata.(strcat('F',num2str(floatnum))) = ...
                rmfield(Mdata.(strcat('F',num2str(floatnum))),mnames{l});
            mnames{l} = [];
            % For informational meta variables
        else
            % Save in metadata structure
            Mdata.(strcat('F',num2str(floatnum))).(names{l}) = ...
                Data.(strcat('F',num2str(floatnum))).(names{l});
            % Remove data fields
            Data.(strcat('F',num2str(floatnum))) = ...
                rmfield(Data.(strcat('F',num2str(floatnum))),names{l});
            names{l} = [];
        end
    end
    % Remove unused variable names
    names = names(~cellfun('isempty',names));
    mnames = mnames(~cellfun('isempty',mnames));
    
    % Add WMO float number to metadata
    Mdata.(strcat('F',num2str(floatnum))).WMO_NUMBER = floatnum;
    
    % CONVERT QUALITY FLAGS TO NUMERIC FORMAT
    for l=1:numel(names)
        if endsWith(names{l},'_QC') && ... % Check for QC identifier
                ~startsWith(names{l},'PROF') % But not a profile QC
            % Vectorize
            Data.(strcat('F',num2str(floatnum))).(names{l}) = ...
                Data.(strcat('F',num2str(floatnum))).(names{l})(:);
            % Replace blanks with zeros
            Data.(strcat('F',num2str(floatnum))).(names{l}) = ...
                strrep(Data.(strcat('F',...
                num2str(floatnum))).(names{l})',' ','0')';
            % Convert to numeric
            Data.(strcat('F',num2str(floatnum))).(names{l}) = ...
                str2num(Data.(strcat('F',num2str(floatnum))).(names{l}));
            % Reshape
            Data.(strcat('F',num2str(floatnum))).(names{l}) = ...
                reshape(Data.(strcat('F',...
                num2str(floatnum))).(names{l}),n_levels,n_prof);
        end
    end
    
    % parse parameter names
    for l=1:numel(mnames)
        if strcmp(mnames{l},'PARAMETER')
            % extract parameter names as coherent strings
            for m = 1:n_param
                temp{m,:} = strrep(Mdata.(strcat('F',...
                    num2str(floatnum))).(mnames{l})(:,m,1,1)',' ','');
            end
            params_keep = ismember(temp,new_vars);
            Mdata.(strcat('F',num2str(floatnum))).(mnames{l}) = ...
                temp(params_keep);
            clear temp;
        end
    end
    
    % parse parameter data modes
    for l=1:numel(mnames)
        if strcmp(mnames{l},'PARAMETER_DATA_MODE')
            % create data mode variable for each parameter
            % expand that variable to match size of data matrix
            z=1;
            for m = 1:n_param
                if params_keep(m)
                    Data.(strcat('F',...
                        num2str(floatnum))).([cell2mat(Mdata.(strcat('F',...
                        num2str(floatnum))).PARAMETER(z)),'_DATA_MODE']) = ...
                        repmat(Mdata.(strcat('F',...
                        num2str(floatnum))).(mnames{l})(m,:),n_levels,1);
                    z=z+1;
                else
                end
            end
        end
    end
    
    % clear both parameter and parameter data mode from metadata
    Mdata.(strcat('F',num2str(floatnum))) = ...
        rmfield(Mdata.(strcat('F',num2str(floatnum))),...
        {'PARAMETER','PARAMETER_DATA_MODE'});
    
    % CONVERT JULD VARIABLE TO SERIAL DATE (SINCE YEAR 1950)
    % AND SAVE AS 'TIME'
    Data.(strcat('F',num2str(floatnum))).('TIME') = ...
        datenum(Data.(strcat('F',num2str(floatnum))).('JULD'))+...
        datenum([1950 1 1]);
    names = [names;'TIME']; % Add 'TIME' to list of variable names
    
    if ~isempty(float_profs)
        for l=1:numel(names)
            % Select only specified profiles
            Data.(strcat('F',...
                num2str(floatnum))).(names{l}) = ...
                Data.(strcat('F',...
                num2str(floatnum))).(names{l})(:,float_profs{n});
        end
    end
end
