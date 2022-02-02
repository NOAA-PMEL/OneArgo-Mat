function success = download_float(floatid, file_type)
% download_float  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   success = download_float(floatid [, file_type])
%
% DESCRIPTION:
%   It downloads the Sprof or meta file for one float with a given floatid.
%
% PREREQUISITE:
%   The Sprof and Meta index files must have been downloaded already.
%
% INPUT:
%   floatid   : WMO ID of a float (integer)
%
% OPTIONAL INPUT:
%   file_type : either 'Sprof' (default) or 'meta'
%
% OUTPUT:
%   success   : 1 for success, 0 for failure
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

global Settings Float Meta;

if nargin < 1
    disp('Usage: download_float(WMO_ID [, file_type])')
    return
end
if nargin < 2
    file_type = 'Sprof';
end

success = 0; % set to 1 after successful download

% make sure Float is initialized
if isempty(Float)
    initialize_argo();
end

if strcmp(file_type, 'Sprof')
    ind = 1:Float.nfloats;
    float_idx = ind(Float.wmoid == floatid);
elseif strcmp(file_type, 'meta')
    ind = 1:length(Meta.wmoid);
    float_idx = ind(Meta.wmoid == floatid);
else
    % so far 'meta' is the only other allowed file type
    warning('unknown file type: %s', file_type)
    return
end

if isempty(float_idx)
    warning('Float %d was not found!', floatid)
    return
end

if strcmp(file_type, 'Sprof')
    local_path = [Settings.prof_dir, Float.file_name{float_idx}];
    url_path = ['dac/', Float.file_path{float_idx}];
    remote_file_update = datenum(Float.update(float_idx), 'yyyymmddHHMMSS');
elseif strcmp(file_type, 'meta')
    local_path = [Settings.meta_dir, Meta.file_name{float_idx}];
    url_path = ['dac/', Meta.file_path{float_idx}];
    remote_file_update = datenum(Meta.update(float_idx), 'yyyymmddHHMMSS');
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
