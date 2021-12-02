function success = download_float(floatid)
% download_float  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   success = download_float(floatid)
%
% DESCRIPTION:
%   It downloads the Sprof file for one float with a given floatid.
%
% PREREQUISITE:
%   The Sprof index file must have been downloaded already. 
%
% INPUT:
%   floatid  : WMO ID of a float (integer)
%
% OUTPUT:
%   success  : 1 for success, 0 for failure
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

global Settings Float;

if nargin < 1
    disp('Usage: download_float(WMO_ID)')
    return
end

success = 0; % set to 1 after successful download

% make sure Float is initialized
if isempty(Float)
    initialize_argo();
end

ind = 1:Float.nfloats;
float_idx = ind(Float.wmoid == floatid);

if isempty(float_idx)
    warning('Float %d was not found!', floatid)
    return
end

local_path = [Settings.prof_dir, Float.file_name{float_idx}];
% now check if the Sprof file exists locally already,
% and if so, if it is up-to-date
if exist(local_path, 'file') == 2
    try
        sprof_date = ncread(local_path, 'DATE_UPDATE')';
        sprof_date = datenum(sprof_date, 'yyyymmddHHMMSS');
        % allow a small tolerance value for numerical imprecision
        if sprof_date > ...
                datenum(Float.update(float_idx), 'yyyymmddHHMMSS') - 0.01
            % existing file has all profiles, no need to download again
            success = 1;
            return;
        end
    catch
        warning('something went wrong, try downloading the file again')
    end
end

success = try_download(['dac/',Float.file_path{float_idx}], ...
    local_path);
