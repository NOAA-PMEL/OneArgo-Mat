function [long_name, units] = get_var_name_units(short_name)
% get_var_name_units  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   [long_name, units] = get_var_name_units(short_name)
%
% DESCRIPTION:
%   This function returns the long name and the units for the variable
%   with the given short name.
%
% INPUT:
%   short_name : case-sensitive name of a variable as it appears in 
%                the Sprof index file, e.g., TEMP or DOXY
%
% OUTPUTS:
%   long_name  : long name of the variable
%   units      : units of the variable
%
% AUTHORS: 
%   H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL), N. Buzby (UW),
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

if nargin < 1
    warning('Usage: get_var_name_units(short_name)')
end

if contains(short_name,'TEMP')
    long_name = 'Temperature';
    units = '(^oC)';
elseif contains(short_name,'PSAL')
    long_name = 'Salinity';
    units = '(PSU)';
elseif contains(short_name,'DENS')
    long_name = 'Density';
    units = '(kg m^{-3})';
elseif contains(short_name,'PRES')
    long_name = 'Pressure';
    units = '(dbar)';
elseif contains(short_name,'DOXY')
    long_name = 'Dissolved Oxygen';
    units = '(\mumol kg^{-1})';
elseif contains(short_name,'NITRATE')
    long_name = 'Nitrate ';
    units = '(\mumol kg^{-1})';
elseif contains(short_name,'AOU')
    long_name = 'Apparent Oxygen Utilization ';
    units = '(\mumol kg^{-1})';
elseif contains(short_name,'PH_IN_SITU')
    long_name = 'pH';
    units = '';
elseif contains(short_name,'CHLA')
    long_name = 'Chlorophyll-a';
    units = '(mg m^{-3})';
elseif contains(short_name,'CDOM')
    long_name = 'Colored Dissolved Organic Matter';
    units = '(ppb)';
elseif contains(short_name,'BBP')
    long_name = 'Backscatter';
    units = '(m^{-1})';
elseif contains(short_name,'DOWN_IRRADIANCE')
    long_name = 'Downwelling irradiance';
    units = '(W m^{-2} nm^{-1})';
elseif contains(short_name,'UP_IRRADIANCE')
    long_name = 'Upwelling irradiance';
    units = '(W m^{-2} nm^{-1})';
elseif contains(short_name,'DOWNWELLING_PAR')
    long_name = 'Downwelling PAR';
    units = '(\mumol Quanta m^{-2} sec^{-1})';
elseif contains(short_name,'BISULFIDE')
    long_name = 'Bisulfide';
    units = '(\mumol kg^{-1})';
elseif contains(short_name,'TURBIDITY')
    long_name = 'Turbidity';
    units = '(ntu)';
else
    warning('unknown variable')
    long_name = [];
    units = [];
end
