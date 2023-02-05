function initialize_sprof(file_name)
% initialize_sprof  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   initialize_sprof(file_name)
%
% DESCRIPTION:
%   This function initializes the global struct Sprof by reading
%   the index file and processing its information.
%
% INPUTS:
%   file_name : name of the index file (with local path)
%
% OUTPUTS: None. Sprof is filled with fields.
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

global Settings Sprof;

fid = fopen(file_name);
H = textscan(fid,'%s %s %f %f %s %d %s %s %s %s','headerlines',9,...
    'delimiter',',','whitespace','');
fclose(fid);
Sprof.urls = H{1};
Sprof.date = H{2};
Sprof.lat  = H{3};
Sprof.lon  = H{4};
Sprof.ocean = H{5};
Sprof.profiler = H{6}; % profiler type
% column 7: institution
Sprof.sens = H{8};
Sprof.split_sens = cellfun(@split, Sprof.sens, 'UniformOutput', false);
Sprof.data_mode = H{9};
Sprof.update = H{10};

% adjust longitude to standard range of -180..180 degrees
Sprof.lon(Sprof.lon > 180) = Sprof.lon(Sprof.lon > 180) - 360;
Sprof.lon(Sprof.lon < -180) = Sprof.lon(Sprof.lon < -180) + 360;

% check for additional (e.g., DOXY2) and unknown sensors
for i = 1:length(Sprof.sens)
    sensors = split(Sprof.sens{i});
    if ~all(ismember(sensors, Settings.avail_vars))
        unknown_sensors = sensors(~ismember(sensors, Settings.avail_vars));
        for s = 1:length(unknown_sensors)
            main_sensor = get_sensor_number(unknown_sensors{s});
            if isempty(main_sensor)
                warning('unknown sensor in index file: %s', unknown_sensors{s})
            else
                Settings.avail_vars{end+1} = unknown_sensors{s};
            end
        end
    end
end

Sprof.wmo = regexp(Sprof.urls,'\d{7}','once','match');
