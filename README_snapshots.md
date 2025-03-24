# Using Argo snapshots

It is recommended to use Argo snapshots for scientific publications. They are released every month and each snapshot comes with a unique DOI. This will ensure that any work based on these data can be reproduced exactly.

However, downloading and unpacking an Argo snapshot is rather time-consuming (several hours for a full snapshot) and takes a lot of disk space, especially for full snapshots. As of January 2024, the tarball for the full snapshot is 63 GB large. Unpacking it with this toolbox takes about 180 GB of disk space temporarily and 70 GB permanently if prof and Sprof nc files are kept for all floats. For BGC snapshots, the disk space requirements are much lower - for the January 2024 snapshot, about 4 GB for the tarball (deleted after unpacking), and 11 GB for the Sprof files.

## Recommended usage

The standard way of using the OneArgo-Mat toolbox relies on files directly downloaded from the GDAC. Only files that are needed by a user will be downloaded, and you will always have the most current versions of them. Therefore, it is recommended to use this method while performing analyses, visualizations etc.

Switching to snapshots is recommended only when the work is being finalized for publication. In order to keep disk usage as low as possible, it is recommended to use the default setting for the 'keep' argument of 0. You can further reduce the amount of required disk space by specifying either the 'dac' or 'floats' arguments to the download_snapshot function. For the 'floats' argument, the return value of the select_profiles function call can be used. (When the 'floats' argument is specified, there is no need to also specify the 'dac' argument, unless you want to further restrict the selection of floats.)

## Notes

### Turn off source control integration

There will be many files in the directory tree of the snapshot during unpacking of full snapshots (mostly individual profile netcdf files, which are not used by this toolbox and will be deleted, unless the 'keep' argument to the download_snapshots function is set to 2 or 3). If the snapshots directory is in the current working directory (or linked to it), the Matlab GUI will by default try to parse these files. This will use up a lot of CPU power and very likely cause an out-of-memory error. To avoid this, turn off source control integration:

In the Matlab GUI under Preferences, under General->Source Control, select "None" instead of "Enable MathWorks source control integration".

### Switching between snapshots and GDAC files

If you want to switch between using snapshot and GDAC files while running Matlab with the toolbox, adjust the value of Settings.use_snapshots in initialize_argo.m and call initialize_argo()

initialize_argo must also be called if the type or date of snapshot files is changed in initialize_argo.m for these changes to be effective.

### Manual downloading and unpacking of snapshot files

You can also download the snapshot file(s) of your choosing directly from
https://www.seanoe.org/data/00311/42182

and untar it yourself from the command line. In case of the full snapshots, you can do this to drastically reduce the amount of required disk space if you only need files from one or few of the DACs. (You can use the list_dacs() function to list the DACs handling the floats of interest.)

Note that the final disk usage is the same if you use the toolbox to unpack the tarball with appropriate options, e.g., for the case shown below:

download_snapshot('dac', 'meds')

Only the temporary disk usage will be higher if you use the download_snapshot command.

However, you must set up the files in a directory structure that is expected by this toolbox.

For instance, if you only need the core float files from the meds DAC (shown here for the September 2023 snapshot):

tar xvzf 104707.tar.gz 202309-ArgoData/dac/meds_core.tar.gz

tar xvzf 104707.tar.gz 202309-ArgoData/dac/ar_index_global_prof.txt.gz

(Optionally, if you don't need any other files from this snapshot: rm 104707.tar.gz)

mkdir Snapshots

mv 202309-ArgoData Snapshots

cd Snapshots/202309-ArgoData/dac

tar xvzf meds_core.tar.gz

mkdir ../Profiles ../Index

mv meds/*/*_prof.nc ../Profiles

gunzip ar_index_global_prof.txt.gz

mv ar_index_global_prof.txt ../Index

cd ..

rm -rf dac





