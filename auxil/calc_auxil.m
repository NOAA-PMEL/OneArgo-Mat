function Data = calc_auxil(Data,varargin)
% calc_auxil  This function is part of the
% MATLAB toolbox for accessing Argo float Data.
%
% USAGE:
%   Data = calc_auxil(Data,varargin)
%
% DESCRIPTION:
%   This function calculates various auxiliary variables from Argo
%   float data: Density, mixed layer depth (MLD), based either on
%   a temperature or a density threshold.
%   Any combination of these variables can be requested.
%
% INPUT:
%   Data      : struct (note that this can be on the original or the
%               interpolated pressure/depth axis)
%
% OPTIONAL INPUTS:
%   'raw',raw                     : use raw data if 'yes', adjusted data
%                                   if no (default: 'no')
%   'calc_dens', calc_dens:         if set, calculate potential density
%   'calc_mld_temp', calc_mld_temp: if set, compute MLD based on T threshold
%   'temp_thresh', temp_threshold : temperature threshold for MLD calculation
%                                   (default: 0.2 dg C); ignored if
%                                   calc_mld_temp is not set to 1
%   'calc_mld_dens', calc_mld_dens: if set, compute MLD based on rho
%                                   threshold
%   'dens_thresh', dens_threshold : density threshold for MLD calculation
%                                   (default: 0.03 kg/m3); ignored if
%                                   calc_mld_dens is not set to 1
%   'float', float_name             float_name (string), e.g.: 'F5904770'
%   (all calc* values are 0=off by default, set them to
%   1=on to activate the calculation)
%
% OUTPUT:
%   Data      : struct with added auxiliary variables - for all variables,
%               either raw or adjusted fields are added to the structure.
%               The raw fields are (adjusted fields have the same names
%               with _ADJUSTED added at the end):
%               calc_dens: DENS (in situ density),
%                          PSAL_ABS (absolute salinity),
%                          TEMP_CNS (conservative temperature)
%               calc_mld_temp: MLD_TEMP (mixed layer depth based on a
%                              temperature threshold)
%               calc_mld_dens: MLD_DENS (mixed layer depth based on a
%                              density threshold)
%
% AUTHORS:
%   J. Sharp, H. Frenzel, A. Fassbender (NOAA-PMEL), N. Buzby (UW)
%
% CITATION:
%   H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2022. OneArgo-Mat:
%   A MATLAB toolbox for accessing and visualizing Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.6588041
%
% LICENSE: oneargo_mat_license.m
%
% DATE: JUNE 1, 2022  (Version 1.0.1)

global Settings;

if nargin < 2
    warning('Usage: calc_auxil(Data,varargin)')
    return
end

% set defaults
raw = 'no';
calc_dens = 0;
calc_mld_temp = 0;
calc_mld_dens = 0;
temp_thresh = Settings.temp_thresh;
dens_thresh = Settings.dens_thresh;
float_name = []; % flag for "unknown"

% parse optional arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'raw')
        raw = varargin{i+1};
    elseif strcmpi(varargin{i}, 'calc_dens')
        calc_dens = varargin{i+1};
    elseif strcmpi(varargin{i}, 'calc_mld_temp')
        calc_mld_temp = varargin{i+1};
    elseif strcmpi(varargin{i}, 'calc_mld_dens')
        calc_mld_dens = varargin{i+1};
    elseif strcmpi(varargin{i}, 'temp_thresh')
        temp_thresh = varargin{i+1};
    elseif strcmpi(varargin{i}, 'dens_thresh')
        dens_thresh = varargin{i+1};
    elseif strcmpi(varargin{i}, 'float_name')
        float_name = varargin{i+1};
    end
end

if ~isempty(float_name)
    float_str = [' for float ', float_name];
else
    float_str = '';
end

% Calculate in situ density:
if calc_dens
    if ~isfield(Data, 'PSAL') && ~isfield(Data, 'PSAL_ADJUSTED')
        warning('No salinity data found, skipping density calculation%s', ...
            float_str);
    else
        try % Try with adjusted values first
            assert(strcmp(raw, 'no'));
            Data.PSAL_ABS_ADJUSTED = gsw_SA_from_SP(Data.PSAL_ADJUSTED,...
                Data.PRES_ADJUSTED,...
                Data.LONGITUDE,Data.LATITUDE); % Calculate absolute salinity
            Data.TEMP_CNS_ADJUSTED = gsw_CT_from_t(Data.PSAL_ABS_ADJUSTED,...
                Data.TEMP_ADJUSTED,...
                Data.PRES_ADJUSTED); % Calculate conservative temperature
            Data.DENS_ADJUSTED = 1000 + gsw_sigma0(Data.PSAL_ABS_ADJUSTED,...
                Data.TEMP_CNS_ADJUSTED); % Calculate potential density
        catch % Then try with unadjusted values
            Data.PSAL_ABS = gsw_SA_from_SP(Data.PSAL,Data.PRES,...
                Data.LONGITUDE,Data.LATITUDE); % Calculate absolute salinity
            Data.TEMP_CNS = gsw_CT_from_t(Data.PSAL_ABS,Data.TEMP,...
                Data.PRES); % Calculate conservative temperature
            Data.DENS     = 1000 + gsw_sigma0(Data.PSAL_ABS,...
                Data.TEMP_CNS); % Calculate in situ density
        end
    end
end

if calc_mld_temp || calc_mld_dens
    try % Try with adjusted values first
        assert(strcmp(raw, 'no'));
        temp = Data.TEMP_ADJUSTED;
    catch % Then try with unadjusted values
        temp = Data.TEMP;
    end
    try % Try with adjusted values first
        assert(strcmp(raw, 'no'));
        pres = Data.PRES_ADJUSTED;
    catch % Then try with unadjusted values
        pres = Data.PRES;
    end
    if ~isfield(Data, 'PSAL') && ~isfield(Data, 'PSAL_ADJUSTED')
        warning('No salinity data found, skipping MLD calculation%s', ...
            float_str);
        calc_mld_temp = 0;
        calc_mld_dens = 0;
    else
        try % Try with adjusted values first
            assert(strcmp(raw, 'no'));
            salt = Data.PSAL_ADJUSTED;
        catch % Then try with unadjusted values
            salt = Data.PSAL;
        end
    end
end

if calc_mld_temp
    % Pre-allocate mixed layer
    Data.MLD_TEMP = nan(1,size(temp,2));
    % Calculate mixed layer depth based on temperature threshold
    for n=1:size(temp,2)
        pressure_prof = pres(:,n); % extract pressure profile
        % determine pressure closest to 10
        [~,ref_idx] = min(abs(pressure_prof-10));
        temperature_prof = temp(:,n); % extract temperature profile
        salt_prof = salt(:,n); % same for salinity
        % compute potential temperature
        ptemp_prof = gsw_pt0_from_t(salt_prof,temperature_prof,pressure_prof);
        % define reference potential temperature as closest to P = 10
        ptemp_ref = ptemp_prof(ref_idx);
        under_ref = pressure_prof > pressure_prof(ref_idx); % index below reference
        % truncate pressure profile to below reference
        pressure_prof = pressure_prof(under_ref);
        % truncate pot. temperature profile to below reference
        ptemp_prof = ptemp_prof(under_ref);
        MLDt_idx = find(ptemp_prof < ptemp_ref-temp_thresh, 1);
        if ~isempty(MLDt_idx)
            Data.MLD_TEMP(1,n) = pressure_prof(MLDt_idx);
        end
    end
end

if calc_mld_dens && (isfield(Data, 'PSAL') || isfield(Data, 'PSAL_ADJUSTED'))
    % Calculate potential density with respect to surface pressure (=0)
    SA = gsw_SA_from_SP(salt,pres,Data.LONGITUDE,Data.LATITUDE);
    CT = gsw_CT_from_t(SA,temp,pres);
    pdensity = 1000 + gsw_sigma0(SA,CT);
    % Pre-allocate mixed layer
    Data.MLD_DENS = nan(1,size(pdensity,2));
    % Identify first instance of potential density that is
    % "dens_thres" greater than surface potential density
    for n = 1:size(pdensity,2)
        pressure_prof = pres(:,n); % extract pressure profile
        % determine pressure closest to 10 dbar; use that as reference
        [~,ref_idx] = min(abs(pressure_prof-10));
        density_prof = pdensity(:,n); % extract n-th density profile
        % define reference density as closest to P = 10 dbar
        density_ref = density_prof(ref_idx);
        if ~isfinite(density_ref)
            % this can happen if the salinity is NaN at 10 dbar
            density_ref = find(isfinite(density_prof(:,1)), 1);
        end
        if isfinite(density_ref) % make sure that a valid density was found
            under_ref = pressure_prof > pressure_prof(ref_idx);
            pressure_prof = pressure_prof(under_ref);
            % truncate pressure profile to below reference
            density_prof = density_prof(under_ref);
            % truncate density profile to below reference
            MLDd_idx = find(density_prof > density_ref+dens_thresh, 1);
            if ~isempty(MLDd_idx)
                Data.MLD_DENS(1,n) = pressure_prof(MLDd_idx);
            end
        end
    end
end
