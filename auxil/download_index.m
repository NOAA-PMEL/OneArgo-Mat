function success = download_index(file_name, file_type)
% download_index  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   success = download_index(file_name, file_type)
%
% DESCRIPTION:
%   This function downloads the specified index file from the GDAC if
%   necessary. If the gzipped file is downloaded, it will be unzipped.
%
% INPUTS:
%   file_name : name of the index file (without path)
%   file_type : type of the index file, e.g. 'Sprof', 'prof', or 'meta'
%
% OUTPUTS:
%   success   : 1 if the file was downloaded successfully or did not
%               need to be downloaded; 0 if it could not be downloaded
%
% AUTHORS:
%   H. Frenzel and J. Sharp (UW-CICOES), A. Fassbender (NOAA-PMEL), N. Buzby (UW)
%
% CITATION:
%   H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2025. OneArgo-Mat:
%   A MATLAB toolbox for accessing and visualizing Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.6588041
%
% LICENSE: oneargo_mat_license.m
%
% DATE: APRIL 16, 2025  (Version 1.1.0)

global Settings;

dest_path = [Settings.index_dir, file_name];

if do_download(dest_path)
    if Settings.verbose
        fprintf('%s index file will now be downloaded.\n', file_type);
    end
    if ~try_download(file_name, dest_path)
        success = 0;
        return;
    end
    if endsWith(file_name, '.gz')
        gunzip(dest_path)
    end
else
    if endsWith(file_name, '.gz')
        gunzip_file = strrep(dest_path, '.gz', '');
        if ~exist(gunzip_file, 'file')
            gunzip(dest_path)
        end
    end
end
success = 1;
