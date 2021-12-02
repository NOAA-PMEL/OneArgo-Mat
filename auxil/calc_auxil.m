function Data = calc_auxil(Data,varargin)
% calc_auxil  This function is part of the
% MATLAB toolbox for accessing BGC Argo float Data.
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
%   'calc_dens', calc_dens:         if set, calculate in situ density
%   'calc_mld_temp', calc_mld_temp: if set, compute MLD based on T threshold
%   'temp_thresh', temp_threshold : temperature threshold for MLD calculation
%                                   (default: 0.2 dg C); ignored if
%                                   calc_mld_temp is not set to 1
%   'calc_mld_dens', calc_mld_dens: if set, compute MLD based on rho
%                                   threshold
%   'dens_thresh', dens_threshold : density threshold for MLD calculation
%                                   (default: 0.03 kg/m3); ignored if
%                                   calc_mld_dens is not set to 1
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
%   J. Sharp, H. Frenzel, A. Fassbender (NOAA-PMEL), N. Buzby (UW),
%   J. Plant, T. Maurer, Y. Takeshita (MBARI), D. Nicholson (WHOI),
%   and A. Gray (UW)
%
% CITATION:
%   H. Frenzel*, J. Sharp*, A. Fassbender, N. Buzby, J. Plant, T. Maurer,
%   Y. Takeshita, D. Nicholson, A. Gray, 2021. BGC-Argo-Mat: A MATLAB
%   toolbox for accessing and visualizing Biogeochemical Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.4971318.
%   (*These authors contributed equally to the code.)
%
% LICENSE: bgc_argo_mat_license.m
%
% DATE: DECEMBER 1, 2021  (Version 1.1)

global Settings;

if nargin < 2
    warning('Usage: calc_auxil(Data,varargin)')
end

% set defaults
calc_dens = 0;
calc_mld_temp = 0;
calc_mld_dens = 0;
temp_thresh = Settings.temp_thresh;
dens_thresh = Settings.dens_thresh;

% parse optional arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'calc_dens')
        calc_dens = varargin{i+1};
    elseif strcmpi(varargin{i}, 'calc_mld_temp')
        calc_mld_temp = varargin{i+1};
    elseif strcmpi(varargin{i}, 'calc_mld_dens')
        calc_mld_dens = varargin{i+1};
    elseif strcmpi(varargin{i}, 'temp_thresh')
        temp_thresh = varargin{i+1};
    elseif strcmpi(varargin{i}, 'dens_thresh')
        dens_thresh = varargin{i+1};
    end
end

% Calculate in situ density:
if calc_dens
    try % Try with adjusted values first
        Data.PSAL_ABS_ADJUSTED = gsw_SA_from_SP(Data.PSAL_ADJUSTED,...
            Data.PRES_ADJUSTED,...
            Data.LONGITUDE,Data.LATITUDE); % Calculate absolute salinity
        Data.TEMP_CNS_ADJUSTED = gsw_CT_from_t(Data.PSAL_ABS_ADJUSTED,...
            Data.TEMP_ADJUSTED,...
            Data.PRES_ADJUSTED); % Calculate conservative temperature
        Data.DENS_ADJUSTED = gsw_rho(Data.PSAL_ABS_ADJUSTED,...
            Data.TEMP_CNS_ADJUSTED,...
            Data.PRES_ADJUSTED); % Calculate in situ density
    catch % Then try with unadjusted values
        Data.PSAL_ABS = gsw_SA_from_SP(Data.PSAL,Data.PRES,...
            Data.LONGITUDE,Data.LATITUDE); % Calculate absolute salinity
        Data.TEMP_CNS = gsw_CT_from_t(Data.PSAL_ABS,Data.TEMP,...
            Data.PRES); % Calculate conservative temperature
        Data.DENS     = gsw_rho(Data.PSAL_ABS,Data.TEMP_CNS,...
            Data.PRES); % Calculate in situ density
    end
end

if calc_mld_temp || calc_mld_dens
    try % Try with adjusted values first
        temp = Data.TEMP_ADJUSTED;
    catch % Then try with unadjusted values
        temp = Data.TEMP;
    end
    try % Try with adjusted values first
        pres = Data.PRES_ADJUSTED;
    catch % Then try with unadjusted values
        pres = Data.PRES;
    end
end
    
if calc_mld_temp
    try % Try with adjusted values first
        sal = Data.PSAL_ADJUSTED;
    catch % Then try with unadjusted values
        sal = Data.PSAL;
    end
    % Pre-allocate mixed layer
    Data.MLD_TEMP = nan(1,size(temp,2));
    % Calculate mixed layer depth based on temperature threshold
    for n=1:size(temp,2)
        pressure_prof = pres(:,n); % extract pressure profile
        % determine pressure closest to 10
        [~,ref_idx] = min(abs(pressure_prof-10));
        temperature_prof = temp(:,n); % extract temperature profile
        sal_prof = sal(:,n); % same for salinity
        % compute potential temperature
        ptemp_prof = gsw_pt0_from_t(sal_prof,temperature_prof,pressure_prof);
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

if calc_mld_dens
    try % Try with adjusted values first
        salt = Data.PSAL_ADJUSTED;
    catch % Then try with unadjusted values
        salt = Data.PSAL;
    end
    % Calculate potential density with respect to surface pressure (=0)
    pdensity = gsw_rho(salt,temp,zeros(size(temp)));
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
