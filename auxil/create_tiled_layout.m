function [ax1, ax2] = create_tiled_layout()
% create_tiled_layout  This function is part of the
% MATLAB toolbox for accessing BGC Argo float data.
%
% USAGE:
%   [ax1, ax2] = create_tiled_layout()
%
% DESCRIPTION:
%   This function creates two tiles with separate axes in the current
%   plotting window.
%
% PREREQUISITE:
%   Plotting window must be open already.
%
% INPUT: None.
%
% OUTPUTS:
%   ax1 : axes for first variable (left and bottom)
%   ax2 : axes for second variable (right and top)
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
% DATE: FEBRUARY 22, 2022  (Version 1.2)

t = tiledlayout(1,1);

ax1 = axes(t);
ax2 = axes(t);

ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax2.XColor = 'b';
ax2.YColor = 'b';

ax1.Box = 'off';
ax2.Box = 'off';

hold(ax2, 'on')
