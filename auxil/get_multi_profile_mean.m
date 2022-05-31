function [mean_prof,std_prof,mean_pres] = get_multi_profile_mean(Datai, ...
    variable)
% get_multi_profile_mean  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   [mean_prof,std_prof,mean_pres] = get_multi_profile_mean(Datai,variable)
%
% DESCRIPTION:
%   This function computes mean and standard deviation of a
%   depth-interpolated variable. Missing values are omitted in the
%   calculations.
%
% INPUTS:
%   Datai     : struct with depth-interpolated fields from multiple floats
%               that must include PRES and the given variable
%   variable  : string with the name of the variable (e.g., DOXY)
%
% OUTPUTS:
%   mean_prof : mean value of the variable across floats (column vector)
%   std_prof  : standard variation of the variable across floats (column vector)
%   mean_pres : mean pressure across floats (column vector)
%
% AUTHORS:
%   H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL), N. Buzby (UW)
%
% CITATION:
%   H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2022. OneArgo-Mat:
%   A MATLAB toolbox for accessing and visualizing Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.6588042
%
% LICENSE: oneargo_mat_license.m
%
% DATE: JUNE 1, 2022  (Version 1.0.1)

floats = fieldnames(Datai);
nfloats = length(floats);

% all Datai structs have the same depth resolution, but the maximum
% value varies by float; need to find the float with the most depths
imax = 0;
max_npres = 0;
this_npres = nan(nfloats,1);
for f = 1:nfloats
    this_npres(f) = size(Datai.(floats{f}).PRES, 1);
    if this_npres(f) > max_npres
        imax = f;
        max_npres = this_npres(f);
    end
end
% this is PRES_ADJUSTED if available:
mean_pres = Datai.(floats{imax}).PRES(:,1); % values are the same in 2nd dim
ndepths = size(mean_pres,1);
total_nprofs = 0;
this_nprofs = nan(nfloats,1);
for f = 1:nfloats
    this_nprofs(f) = size(Datai.(floats{f}).PRES,2);
    total_nprofs = total_nprofs + this_nprofs(f);
end

all_profs = nan(ndepths, total_nprofs); % pre-allocate
count_profs = 0;
for f = 1:nfloats
    all_profs(1:this_npres(f),...
        count_profs+1:count_profs + this_nprofs(f)) = ...
        Datai.(floats{f}).(variable);
    count_profs =  count_profs + this_nprofs(f);
end
mean_prof = mean(all_profs,2,'omitnan');
std_prof = std(all_profs,[],2,'omitnan');
