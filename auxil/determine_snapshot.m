function determine_snapshot()
% determine_snapshot  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   determine_snapshot()
%
% DESCRIPTION:
%   This function determines which snapshot of Argo data will be used,
%   based on the values of Settings.use_snapshots and
%   Settings.default_type.
%   Settings.snap_path, Settings.snap_date, Settings.snap_url, 
%   Settings.snap_file, and Settings.snap_size will be modified.
%
% INPUT: None
%
% OUTPUT: None.
%   Global variable Settings will be modified.     
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

global Settings;

snap_date = Settings.use_snapshots; % shorthand names
snap_type = Settings.default_type;
Settings.snap_path = [];

if ~snap_date
    return
elseif snap_date > 1 && (snap_date < 201212 || snap_date > 203800)
    warning('Wrong format of snapshot date, must be 1 or YYYYMM')
    return
end

% set search patterns for parsing of web page
if strcmp(snap_type, 'all') || strcmp(snap_type, 'phys')
    pattern = 'Global GDAC Argo data files (';
    pat_month = 'Global GDAC Argo data files \(20\d{2}-\d{2}-\d{2} snapshot';
elseif strcmp(snap_type, 'bgc')
    pattern = 'BGC Sprof data files (';
    pat_month = 'BGC Sprof data files \(20\d{2}-\d{2}-\d{2} snapshot\)';
else
    fprintf('No such type of snapshot file: "%s"\n', snap_type);
    return
end

mth_idx = length(pattern) + 1; % first character of YYYY-MM

% download the web page that contains the links to all snapshots
try
    page = webread('https://www.seanoe.org/data/00311/42182/');
    idx = strfind(page, pattern);
catch
    if check_existing_snapshot()
        disp('The page at seanoe.org could not be reached,')
        disp('but the requested snapshot appears to exist already in:')
        disp(Settings.snap_path)
        return;
    else
        idx = [];
    end
end
if isempty(idx)
    fprintf('The snapshot web page could not be read properly.\n')
    return
end

% determine which snapshots are available
nidx = length(idx);
snap_month = nan(nidx, 1);
snap_size = nan(nidx, 1);
snap_url = cell(nidx, 1);
count = 0;
for i = 1:nidx
    this_link = page(idx(i):idx(i)+400);
    match_month = regexp(this_link, pat_month, 'match', 'once');
    pat_url = '"fileUrl":"http[\w/:.]+"';
    match_url = regexp(this_link, pat_url, 'match', 'once');
    pat_size = '"size":\d+';
    match_size = regexp(this_link, pat_size, 'match', 'once');
    % note that some snapshots are only available on demand, for those,
    % there is no fileUrl entry
    if ~isempty(match_month) && ~isempty(match_url)
        count = count + 1;
        year = str2double(match_month(mth_idx:mth_idx+3));
        month = str2double(match_month(mth_idx+5:mth_idx+6));
        snap_month(count) = 100 * year + month; % YYYYMM
        snap_url{count} = match_url(12:end-1);
        if ~isempty(match_size)
            snap_size(count) = uint64(str2double(match_size(8:end)));
        else
            snap_size(count) = -1;
        end
    end
end

if ~count
    fprintf('No matching snapshots were found.\n')
    return
end

if snap_date == 1
    % use the most recent snapshot, need to sort the available ones by snap_date
    % most recent will be first
    [~,isort] = sort(snap_month(isfinite(snap_month)), 'descend');
    isnap = isort(1);
    Settings.snap_date = snap_month(isnap);
else
    % check if the specified snapshot exists
    [mini,isnap] = min(abs(snap_month - snap_date));
    if mini > 0
        fprintf('No snapshot found for %d\n', snap_date)
        fprintf('The nearest available snapshot is %d\n', snap_month(isnap))
        return
    end
    Settings.snap_date = snap_date;
end

if strcmp(snap_type, 'all') || strcmp(snap_type, 'phys')
    Settings.snap_path = sprintf('%d-ArgoData/', Settings.snap_date);
else
    Settings.snap_path = sprintf('%d-BgcArgoSprof/', Settings.snap_date);
end

Settings.snap_url = snap_url{isnap};
% this is the most commonly used format for file names:
Settings.snap_file = regexp(snap_url{isnap}, '\d+\.tar.gz', 'match', 'once');
if isempty(Settings.snap_file) % alternate file name format
    Settings.snap_file = regexp(snap_url{isnap}, '\d+\.tgz', 'match', 'once');
end
if isempty(Settings.snap_file) % this should not happen
    fprintf('unexpected file name for snapshot: "%s"', snap_url{isnap})
    return
end

Settings.snap_size = snap_size(isnap); % in GB
