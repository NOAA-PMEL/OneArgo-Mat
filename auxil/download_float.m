function success = download_float(floatid, file_type)
% download_float  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   success = download_float(floatid [, file_type])
%
% DESCRIPTION:
%   It downloads the Sprof or meta file for one float with a given floatid.
%
% PREREQUISITE:
%   The Sprof, Prof, Meta, and Tech index files must have been
%   downloaded already.
%
% INPUT:
%   floatid   : WMO ID of a float (integer)
%
% OPTIONAL INPUT:
%   file_type : 'prof' (default), 'Sprof', 'meta', or 'tech'
%
% OUTPUT:
%   success   : 1 for success, 0 for failure
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

global Settings Float Meta Tech Traj;

if nargin < 1
    disp('Usage: download_float(WMO_ID [, file_type])')
    return
end
if nargin < 2
    file_type = 'prof';
end

success = 0; % set to 1 after successful download

% make sure Float is initialized
if isempty(Float)
    initialize_argo();
end

if contains(file_type, 'prof')
    ind = 1:Float.nfloats;
    float_idx = ind(Float.wmoid == floatid);
elseif strcmp(file_type, 'meta')
    ind = 1:length(Meta.wmoid);
    float_idx = ind(Meta.wmoid == floatid);
elseif strcmp(file_type, 'tech')
    ind = 1:length(Tech.wmoid);
    float_idx = ind(Tech.wmoid == floatid);
elseif strcmp(file_type, 'traj')
    ind = 1:length(Traj.wmoid);
    float_idx = ind(Traj.wmoid == floatid);
    if length(float_idx) > 1
        float_idx = float_idx(contains(Traj.file_name(float_idx),'Dtraj'));
    end
else
    warning('unknown file type: %s', file_type)
    return
end

if isempty(float_idx)
    warning('Float %d was not found!', floatid)
    return
end

if contains(file_type, 'prof')
    local_path = [Settings.prof_dir, Float.file_name{float_idx}];
    url_path = ['dac/', Float.file_path{float_idx}];
    remote_file_update = datenum(Float.update(float_idx), 'yyyymmddHHMMSS');
elseif strcmp(file_type, 'meta')
    local_path = [Settings.meta_dir, Meta.file_name{float_idx}];
    url_path = ['dac/', Meta.file_path{float_idx}];
    remote_file_update = datenum(Meta.update(float_idx), 'yyyymmddHHMMSS');
elseif strcmp(file_type, 'tech')
    local_path = [Settings.tech_dir, Tech.file_name{float_idx}];
    url_path = ['dac/', Tech.file_path{float_idx}];
    remote_file_update = datenum(Tech.update(float_idx), 'yyyymmddHHMMSS');
elseif strcmp(file_type, 'traj')
    local_path = [Settings.traj_dir, Traj.file_name{float_idx}];
    url_path = ['dac/', Traj.file_path{float_idx}];
    remote_file_update = datenum(Traj.update(float_idx), 'yyyymmddHHMMSS');
end

% now check if the specified file exists locally already,
% and if so, if it is up-to-date
if exist(local_path, 'file') == 2
    try
        local_file_update = ncread(local_path, 'DATE_UPDATE')';
        local_file_update = datenum(local_file_update, 'yyyymmddHHMMSS');
        % allow a small tolerance value for numerical imprecision
        if local_file_update > remote_file_update - 0.1
            % existing file is up-to-date, no need to download again
            success = 1;
            return;
        end
    catch
        warning('something went wrong, try downloading the file again')
    end
end

success = try_download(url_path, local_path);
