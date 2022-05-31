function Datai = depth_interp(Data, qc_flags, varargin)
% depth_interp  This function is part of the
% MATLAB toolbox for accessing BGC Argo float Data.
%
% USAGE:
%   Datai = depth_interp(Data, qc_flags, varargin)
%
% DESCRIPTION:
%   This function interpolates values for BGC-Argo parameters against depth
%   for a set of float profiles.
%
% INPUTS:
%   Data     : struct with data for one or more floats; expected fields are
%              PRES, TIME, LONGITUDE, LATITUDE, CYCLE_NUMBER
%   qc_flags : array with accepted QC flags (see plot_profiles.m for a full
%              description)
%
% OPTIONAL INPUTS:
%   'prs_res',prs_res             : pressure resolution (default: 2 dbar)
%   'raw',raw                     : use raw data if 'yes', adjusted data
%                                   if no (default: 'no')
%   'calc_dens',calc_dens         : if set to 1, calculate potential density on
%                                   interpolated depth levels
%   'calc_mld_dens',calc_mld_dens : if set to 1, calculate mixed layer
%                                   depth (MLD) based on a density criterion
%   'dens_thres',dens_thres       : density threshold for MLD calculation;
%                                   default value is set in initialize_argo
%   'calc_mld_temp',calc_mld_temp : if set to 1, calculate mixed layer
%                                   depth based on a temperature criterion
%   'dens_thres',dens_thres       : temperature threshold for MLD calculation;
%                                   default value is set in initialize_argo
%   (Note that MLD can be computed both ways at the same time.)
%
% OUTPUT:
%   Datai : struct with depth-interpolated variables
%
% AUTHORS:
%   J. Sharp, H. Frenzel, A. Fassbender (NOAA-PMEL), N. Buzby (UW)
%
% CITATION:
%   H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2022. OneArgo-Mat:
%   A MATLAB toolbox for accessing and visualizing Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.6588042
%
% LICENSE: oneargo_mat_license.m
%
% DATE: JUNE 1, 2022  (Version 1.0.1)

global Settings;

if nargin < 2
    warning('Usage: depth_interp(Data, qc_flags, [, varargin])')
    return
end

% set defaults
prs_res = 2;
calc_dens = 0;
calc_mld_dens = 0;
calc_mld_temp = 0;
varargpass= {};
raw = 'no';

% parse optional arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'prs_res')
        prs_res = varargin{i+1};
    elseif strcmpi(varargin{i}, 'raw')
        raw = varargin{i+1};
    elseif strcmpi(varargin{i}, 'calc_dens')
        calc_dens = varargin{i+1};
    elseif strcmpi(varargin{i}, 'calc_mld_dens')
        calc_mld_dens = varargin{i+1};
    elseif strcmpi(varargin{i}, 'calc_mld_temp')
        calc_mld_temp = varargin{i+1};
    elseif strcmpi(varargin{i}, 'temp_thresh') || ...
            strcmpi(varargin{i}, 'dens_thresh')
        varargpass = [varargpass, varargin{i:i+1}];
    end
end

% DEFINE PRESSURE DATA AS 'X'
try
    assert(strcmp(raw, 'no'));
    X = Data.PRES_ADJUSTED;
    good_vals = sum(isfinite(X));
    % a somewhat arbitrary criterion: at least half of the profiles
    % must have valid adjusted pressure values, or switch to
    % raw pressure
    assert(sum(good_vals) > 0.5 * length(good_vals));
catch
    X = Data.PRES;
end

% CONSTRUCT DEPTH AXIS ON WHICH TO INTERPOLATE DEPENDENT VARIABLES
xi = (0:prs_res:prs_res*ceil(max(max(X))/prs_res))';
xi = repmat(xi,1,size(X,2));

% START LOOP FOR EACH DEPENDENT VARIABLE
vars = fieldnames(Data);
for k=1:numel(vars)
    % CHECK FOR EXISTENCE OF FIELD
    if startsWith(vars{k},Settings.avail_vars) && ...
            ~startsWith(vars{k}, 'PRES') && ~endsWith(vars{k},'_QC')
        % DEFINE DEPENDENT VARIABLE AS 'Y'
        Y = Data.(vars{k});
        
        % PRE-ALLOCATE INTERPOLATED DEPENDENT VARIABLE
        yi = NaN(size(xi,1),size(xi,2));
        
        % LOOP INTERPOLATIONS FOR EACH PROFILE
        for n = 1:size(X,2)
            if ~sum(isfinite(X(:,n))) || ~sum(isfinite(Y(:,n)))
                yi(:,n)=NaN; % FILL VECTOR WITH NaN IF EMPTY INPUT
            else
                % INTERPOLATE FOR EACH PROFILE
                if isfield(Data, [vars{k}, '_QC'])
                    qc_mask = ismember(Data.([vars{k},'_QC'])(:,n),...
                        qc_flags);
                else
                    qc_mask = ones(size(X(:,1)));
                end
                idx = ~isnan(X(:,n)) & ~isnan(Y(:,n)) & qc_mask;
                x_filt  = X(:,n);
                x_filt = x_filt(idx);
                y_filt = Y(:,n);
                y_filt = y_filt(idx);
                try
                    yi(:,n) = interp1(x_filt,y_filt,xi(:,n));
                catch
                end
            end
        end
        Datai.(vars{k}) = yi; % ADD INTERPOLATED DEPENDENT VARIABLE TO OUTPUT
    end
end

Datai.PRES = xi; % ADD INTERPOLATED PRESSURE GRID TO OUTPUT

% RESHAPE AND ADD DATE, LATITUDE, AND LONGITUDE TO OUTPUT
Datai.TIME = repmat(Data.TIME(1,:),size(xi,1),1);
Datai.LATITUDE = repmat(Data.LATITUDE(1,:),size(xi,1),1);
Datai.LONGITUDE = repmat(Data.LONGITUDE(1,:),size(xi,1),1);
Datai.CYCLE_NUMBER = repmat(Data.CYCLE_NUMBER(1,:),size(xi,1),1);

if calc_dens || calc_mld_dens || calc_mld_temp
    Datai = calc_auxil(Datai, 'calc_dens', 1, ...
        'calc_mld_dens', calc_mld_dens, 'calc_mld_temp', calc_mld_temp, ...
        'raw', raw, varargpass{:});
end
