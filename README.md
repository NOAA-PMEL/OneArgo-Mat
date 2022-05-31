# Argo Toolbox for MATLAB (OneArgo-Mat)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6588041.svg)](https://doi.org/10.5281/zenodo.6588041)


## NOTE

This is the first release of this toolbox under the OneArgo-Mat name. It is a further development of the [BGC-Argo-Mat](https://github.com/NOAA-PMEL/BGC_Argo_Mat_Toolbox) toolbox.

## ABOUT

This toolbox contains a variety of functions for accessing, processing, and visualizing [Argo](https://argo.ucsd.edu) data, including Core, Deep, and Biogeochemical Argo. Functions are designed to be maximally efficient, to provide access to the most up-to-date data available, and to allow for downloading and plotting of those data based on numerous user-defined conditions.

## QUICK DEMO

This <a href="Demo.md">sample code</a> shows an example of selecting BGC floats that match geographic and temporal limits and visualizing some of the data.

## INSTALLATION AND USE

This repository can be cloned using the command line or the GitHub Desktop application. The files can also be downloaded directly in zipped format.

Before use, make sure the files from this repository are placed in a directory that is in the MATLAB search path. Or add the directory where they are located to the search path (https://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html). Or run it in the directory where the Tutorial.m script was placed.

For an overview of how to use this toolbox, step through the Tutorial script (using the "Run and Advance" button of the MATLAB editor window), a tutorial that was developed for the [2021 GO-BGC Scientific Workshop](https://www.us-ocb.org/joint-gobgc-workshop/) (6/30/21).

## FUNCTIONS

### Main functions (to be called from script or command window):

download_meta_files.m    : downloads <WMOID>_meta.nc files from GDAC<br/>
download_tech_files.m    : downloads <WMOID>_tech.nc files from GDAC<br/>
download_traj_files.m    : downloads <WMOID>_traj.nc files from GDAC<br/>
get_lon_lat_time.m       : returns longitude, latitude, time for selected floats and profiles<br/>
initialize_argo.m        : defines standard settings and paths and downloads synthetic profile index file<br/>
list_sensors.m           : shows available sensors across fleet or for specified floats<br/>
load_float_data.m        : loads data of one or more specified float(s) into memory<br/>
qc_filter.m              : filters variables in a Data structure by QC flags<br/>
select_profiles.m        : returns profiles and corresponding floats based on input criteria<br/>
show_profiles.m          : downloads float data and calls plot_profiles to create plot<br/>
show_sections.m          : downloads float data and calls plot_sections to create plot<br/>
show_timeseries.m        : downloads float data and calls plot_timeseries to create plot<br/>
show_trajectories.m      : downloads float data and calls plot_trajectories to create plot<br/>

Use "help function_name" in the MATLAB command window to see a full description of input and output parameters for these functions.

### Background functions (primarily called by main functions in background):

calc_auxil.m             : calculates various auxiliary variables from Argo float data<br/>
check_datenum.m          : checks if input matches datenum format<br/>
check_dir.m              : determines if a directory needs to be created and does so if necessary<br/>
check_float_variables.m  : checks if specified variables are available for the specified floats<br/>
check_variables.m        : checks if specified variables are available<br/>
combine_variables.m      : creates a cell array with variables names (basic and extended sets)<br/>
create_tiled_layout.m    : creates tiled layout in a profile plot with two different variables<br/>
depth_interp.m           : interpolates values for BGC-Argo parameters against depth<br/>
do_download.m            : determines if a file should be downloaded or not<br/>
do_pause.m               : pauses execution of main_workshop (if used without desktop)<br/>
download_float.m         : downloads the prof, Sprof or other NetCDF file for one float from the GDAC<br/>
download_index.m         : downloads one index file from the GDAC<br/>
download_multi_floats.m  : calls download_float to download NetCDF files from the GDAC for multiple floats<br/>
get_dims.m               : returns number of profiles, parameters, and depth levels of one Sprof file<br/>
get_inpolygon.m          : determines which of the given points are within the specified lon/lat limits<br/>
get_lon_lat_lims.m       : obtains maximum/minimum latitude and longitude values from input data<br/>
get_multi_profile_mean.m : calculates the mean profile of multiple profiles<br/>
get_sensor_number.m      : determines sensor name and number from input<br/>
get_var_name_units.m     : returns the long variable name and units name for a given short parameter name input<br/>
initialize_meta.m        : reads the meta index file and initializes the global Meta struct<br/>
initialize_prof.m        : reads the prof index file and initializes the global Prof struct<br/>
initialize_sprof.m       : reads the sprof index file and initializes the global Sprof struct<br/>
initialize_tech.m        : reads the tech index file and initializes the global Tech struct<br/>
initialize_traj.m        : reads the traj index file and initializes the global Traj struct<br/>
interp_lonlat.m          : interpolates missing positions if possible<br/>
mod_xaxis_time.m         : sets length and labels for x axis (time) in time series and section plots<br/>
plot_one_profile.m       : plots one profile<br/>
plot_profiles.m          : plots profiles of one or more specified float(s) for the specified variable(s)<br/>
plot_sections.m          : plots sections of one or more specified float(s) for the specified variable(s)<br/>
plot_timeseries.m        : plots time series of one or more specified float(s)<br/>
plot_trajectories.m      : plots trajectories of one or more specified float(s)<br/>
select_profiles_per_type.m : helper function of select_profiles that selects one type of floats<br/>
set_xlim.m               : sets xlimits of current plot (if start/end dates were specified)<br/>
try_download.m           : attempts to download a file from any of the specified GDACs<br/>


## DRIVER SCRIPTS
Tutorial.m               : tutorial script for GO-BGC Scientific Workshop (6/30/2021)<br/>
Comp_doxy_doxy2.m        : shows use of multiple sensors in select_profiles, show_timeseries, show_profiles<br/>
Demo.mlx                 : live script version of the third example from Tutorial.m<br/>
Show_deep_floats.m       : shows use of 'depth' option for select_profiles, calls list_sensors<br/>
Show_doxy_rmode.m        : shows use of 'mode' and 'min_num_prof' options for select_profiles<br/>
Show_missing_lonlat.m    : shows interpolation of missing locations in select_profiles and show_trajectories<br/>
Show_under_ice.m         : shows marking of estimated positions in show_trajectories<br/>

## REQUIREMENTS
At least MATLAB R2016b (or any newer release) is needed to use these functions without modifications.<br/>
At least MATLAB R2019b is needed to overlay profiles of two different types of variables.<br/>
An Internet connection is needed to get the latest versions of index and profile files; but analysis and visualization functions can be run offline.<br/>
Memory requirements depend on the number of profiles and variables that are simultaneously loaded into memory.<br/>
Up to 10 GB local disk space is needed to store all Sprof files of BGC floats.<br/>
Up to 50 GB local disk space is needed to store all prof files of core and deep floats.<br/>
<br/>
Downloading all Sprof files of BGC floats for the first time takes about 30 minutes with a fast Internet connection.<br/>
Downloading all prof files of core and deep floats for the first time will take at least four hours, even with a fast Internet connection.


## COMMENTS, BUGS etc.?
Please feel free to use the GitHub Issues and Pull Requests features to report any problems with this code and to suggest bug fixes or additional features.

## BGC-ARGO GUIDE
More detailed information about quality control flags, raw and adjusted modes, etc., can be found in
H. C. Bittig et al., 2019. Front. Mar. Sci. https://doi.org/10.3389/fmars.2019.00502.

## TOOLBOX IN OTHER LANGUAGES
This toolbox has been translated to R:<br/>
[R toolbox](https://github.com/NOAA-PMEL/OneArgo-R)

A similar toolbox in Python:<br/>
[Python toolbox](https://github.com/go-bgc/workshop-python)

[Video tutorials for the toolbox](https://www.go-bgc.org/getting-started-with-go-bgc-data)

## CITATION

Please cite this toolbox as:

H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2022. OneArgo-Mat: A MATLAB toolbox for accessing and visualizing Argo data. Zenodo. https://doi.org/10.5281/zenodo.6588041

## LEGAL DISCLAIMER

This repository is a software product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA), or the United States Department of Commerce (DOC). All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the DOC or DOC bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the DOC. The DOC seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by the DOC or the United States Government.
