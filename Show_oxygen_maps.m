%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show_oxygen_maps.m
%
% This example script for the OneArgo-Mat toolbox finds floats with oxygen
% data in the Eastern Pacific that have undergone delayed mode quality
% control (DMQC) and demonstrates the plotting of maps at specified depths.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lon_lim = [-150 -65]; 
lat_lim = [-40 40]; 
start_date = [2018,1,1]; 
end_date = [];

% select floats and profiles based on criteria
[doxy_floats,doxy_float_profs] = select_profiles(lon_lim, lat_lim, ...
    start_date, end_date, ...
    'sensor', 'DOXY', 'mode', 'D', 'ocean', 'P');
disp(['# of matching profiles: ' num2str(sum(cellfun('length',...
    doxy_float_profs)))]);
disp(['# of matching floats: ' num2str(length(doxy_floats))]);

% show all trajectories
show_trajectories(doxy_floats,'color','multiple',...
    'float_profs',doxy_float_profs,'legend','no','png','DOXY_Dmode_EPac.png');

% show oxygen maps at 20 and 800 dbar depths
show_maps(doxy_floats,{'DOXY'},[20,800],...
    'float_profs',doxy_float_profs,'png','Map','caxis',[0 300],'raw','no_strict');

% load data for further analysis
% Data = load_float_data(doxy_floats,{'TEMP';'DOXY'},doxy_float_profs);




