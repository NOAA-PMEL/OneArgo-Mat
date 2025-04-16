function mod_xaxis_time(min_time, max_time, start_date, end_date, ...
    time_label)
% mod_xaxis_timeseries  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   mod_xaxis_time(min_time, max_time, start_date, end_date, time_label)
%
% DESCRIPTION:
%   This function modifies the x axis of a time series or section plot by
%   setting its limits first to [min_time, max_time] and then applies
%   non-empty values of start_date and end_date. If time_label is empty,
%   the length of the plot determines the choice of time label -
%   years for plots of at least 1.5 years, months for plots of over 60 days,
%   days otherwise.
%
% PREREQUISITE:
%   The time series or section plot must exist already.
%
% INPUTS:
%   min_time   : datenum value for the first time value of the data
%   max_time   : datenum value for the last time value of the data
%   start_date : datenum value for the starting point (may be empty)
%   end_date   : datenum value for the ending point (may be empty)
%   time_label : type of time label: either years ('y'), months ('m'),
%                or days ('d'); if left empty, default is used,
%                which depends on length of time shown:
%                'd' for up to 60 days, 'm' for up to 18 months,
%                'y' otherwise
%
% OUTPUT: None
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

if max_time > min_time
    xlim([min_time, max_time]); % tight layout
end
xl = set_xlim(start_date, end_date);

% determine type of time label based on length of time series
if isempty(time_label)
    if xl(2) - xl(1) >  548 % in days; ~1.5 years
        time_label = 'y';
    elseif xl(2) - xl(1) > 60
        time_label = 'm';
    else
        time_label = 'd';
    end
end
if strncmpi(time_label, 'y', 1)
    set(gca,'XTick',datenum([(2000:2030)' ones(31,1) ones(31,1)]));
    datetick('x','yyyy','keeplimits','keepticks');
    xlabel('Year','FontSize',14);
elseif strncmpi(time_label, 'm', 1)
    set(gca,'XTick',datenum([repelem((2000:2030)',12,1) ...
        repmat((1:12)',31, 1) ones(31*12,1)]));
    datetick('x','mm-yyyy','keeplimits','keepticks');
    xlabel('Month','FontSize',14);
else
    set(gca,'XTick',datenum(2000,1,1):3:datenum(2030,12,31));
    datetick('x','mm-dd-yyyy','keeplimits','keepticks');
    xlabel('Day','FontSize',14);
end
