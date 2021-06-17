function success = try_download(filename, dest_path)
% try_download  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   try_download(filename, dest_path)
%
% DESCRIPTION:
%   This function attempts to download a file from any of the GDACs
%   specified in the Settings.hosts cell array.
%
% INPUTS:
%   filename  : name of the file at the GDAC
%   dest_path : full (relative or absolute) path to the local file
%
% OUTPUT:
%   success   : 1 for successul download; 2 for unsuccessful download,
%               but the file exists locally already; 0 for failure
%
% AUTHORS: 
%   H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL),
%   J. Plant, T. Maurer, Y. Takeshita (MBARI), D. Nicholson (WHOI),
%   and A. Gray (UW)
%
% CITATION:
%   H. Frenzel*, J. Sharp*, A. Fassbender, J. Plant, T. Maurer,
%   Y. Takeshita, D. Nicholson, A. Gray, 2021. BGC-Argo-Mat: A MATLAB
%   toolbox for accessing and visualizing Biogeochemical Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.4971318.
%   (*These authors contributed equally to the code.)
%
% LICENSE: bgc_argo_mat_license.m
%
% DATE: June 15, 2021

global Settings;

success = 0; % default: failure
for h = 1:length(Settings.hosts)
    try
        if Settings.verbose
            fprintf('Attempting download of %s\nto %s ... ', ...
                [Settings.hosts{h}, filename], dest_path);
        end
        websave(dest_path, [Settings.hosts{h}, filename]);
        if Settings.verbose
            fprintf('success!\n');
        end
        success = 1;
        break
    catch
        if Settings.verbose
            fprintf('failure!\n');
        end
        if h == length(Settings.hosts)
            if exist(dest_path, 'file') ~= 2
                return
            else
                % caller may distinguish between return values 1 and 2
                success = 2; 
            end
        end
    end
end
