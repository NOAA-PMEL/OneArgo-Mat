# BGC-Argo Toolbox for MATLAB (BGC-Argo-Mat)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4971318.svg)](https://doi.org/10.5281/zenodo.4971318)

## ABOUT

This toolbox contains a variety of functions for accessing, processing, and visualizing [Biogeochemical Argo](https://biogeochemical-argo.org) data. Functions are designed to be maximally efficient, to provide access to the most up-to-date data available, and to allow for downloading and plotting of those data based on numerous user-defined conditions.

## INSTALLATION AND USE

This repository can be cloned using the command line or the GitHub Desktop application. The files can also be downloaded directly in zipped format.

Before use, make sure the files from this repository are placed in a directory that is in the MATLAB search path. Or add the directory where they are located to the search path (https://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html). Or run it in the directory where the main_workshop script was placed.

For an overview of how to use this toolbox, step through the 'main_workshop' script (using the "Run and Advance" button of the MATLAB editor window), a tutorial that was developed for the [2021 GO-BGC Scientific Workshop](https://www.us-ocb.org/joint-gobgc-workshop/) (6/30/21).

## FUNCTIONS

### Main functions (to be called from script or command window):

get_lon_lat_time.m       : returns longitude, latitude, time for selected floats and profiles<br/>
initialize_argo.m        : defines standard settings and paths and downloads synthetic profile index file<br/>
load_float_data.m        : loads data of one or more specified float(s) into memory<br/>
qc_filter.m              : filters variables in a Data structure by QC flags<br/>
select_profiles.m        : returns profiles and corresponding floats based on input criteria<br/>
show_profiles.m          : downloads float data and calls plot_profiles to create plot<br/>
show_sections.m          : downloads float data and calls plot_sections to create plot<br/>
show_trajectories.m      : downloads float data and calls plot_trajectories to create plot<br/>
main_workshop.m          : tutorial script for GO-BGC Scientific Workshop (6/30/21)<br/>

### Background functions (primarily called by main functions in background):
calc_auxil.m             : calculates various auxiliary variables from Argo float data<br/>
check_dir.m              : determines if a directory needs to be created and does so if necessary<br/>
combine_variables.m      : creates a cell array with variables names (basic and extended sets)<br/>
depth_interp.m           : interpolates values for BGC-Argo parameters against depth<br/>
do_download.m            : determines if a file should be downloaded or not<br/>
do_pause.m               : pauses execution of main_workshop (if used without desktop)<br/>
download_float.m         : downloads the Sprof NetCDF file for one float<br/>
download_multi_floats.m  : calls download_float to download Sprof NetCDF files for multiple floats<br/>
get_dims.m               : returns number of profiles, parameters, and depth levels of one Sprof file<br/>
get_inpolygon.m          : determines which of the given points are within the specified lon/lat limits<br/>
get_lon_lat_lims.m       : obtains maximum/minimum latitude and longitude values from input data<br/>
get_multi_profile_mean   : calculates the mean profile of multiple profiles<br/>
get_var_name_units.m     : returns the long variable name and units name for a given short parameter name input<br/>
plot_profiles.m          : plots profiles of one or more specified float(s) for the specified variable(s)<br/>
plot_sections.m          : plots sections of one or more specified float(s) for the specified variable(s)<br/>
plot_trajectories.m      : plots trajectories of one or more specified float(s)<br/>
try_download.m           : attempts to download a file from any of the specified GDACs<br/>

Use "help function_name" in the MATLAB command window to see a full description of input and output parameters for these functions.

## REQUIREMENTS
At least MATLAB R2016b (or any newer release) is needed to use these functions without modifications.<br/>
An Internet connection is needed to get the latest versions of index and Sprof files; but the repository includes versions of these files so that it can be run offline.<br/>
Memory requirements depend on the number of profiles and variables that are simultaneously loaded into memory.

## COMMENTS, BUGS etc.?
Please feel free to use the GitHub Issues and Pull Requests features to report any problems with this code and to suggest bug fixes or additional features.

## BGC-ARGO GUIDE
More detailed information about quality control flags, raw and adjusted modes, etc., can be found in
H. C. Bittig et al., 2019. Front. Mar. Sci. https://doi.org/10.3389/fmars.2019.00502.


## TOOLBOX IN OTHER LANGUAGES
This toolbox has been translated to R:<br/>
[R toolbox](https://github.com/euroargodev/BGC-ARGO_R_WORKSHOP)

A similar toolbox in Python:<br/>
[Python toolbox](https://github.com/go-bgc/workshop-python)

[Video tutorials for the toolbox](https://www.go-bgc.org/getting-started-with-go-bgc-data)

## CITATION

Please cite this toolbox as:

H. Frenzel*, J. Sharp*, A. Fassbender, N. Buzby, J. Plant, T. Maurer, Y. Takeshita, D. Nicholson, A. Gray, 2021. BGC-Argo-Mat: A MATLAB toolbox for accessing and visualizing Biogeochemical Argo data. Zenodo. https://doi.org/10.5281/zenodo.4971318.

*These authors contributed equally to the code.

## LEGAL DISCLAIMER

This repository is a software product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA), or the United States Department of Commerce (DOC). All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the DOC or DOC bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the DOC. The DOC seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by the DOC or the United States Government.
