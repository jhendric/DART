% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

% Select field to plot (U, V, W, GZ, T, MU, QV, QC, QR)

     field_num = input('Input field type, 1=U, 2=V, 3=W, 4=GZ, 5=T, 6=MU, 7=QV, 8=QC, 9=QR: ');

% Get file name of true state file
fname = 'True_State';
tlon = getnc(fname, 'XLON');
we = size(tlon, 2);
tlat = getnc(fname, 'XLAT');
sn = size(tlat, 1);
level = getnc(fname, 'level');
bt = size(level, 1);

state_vec = getnc(fname, 'state');

% Get a time level from the user
itime = input('Input time level: ');

single_state = state_vec(itime, :);

% Get level for free atmosphere fields
if field_num == 6
   field_level = 1;
else
   field_level = input('Input level: ');
end

start_var = 1
nx = we + 1
ny = sn
var_units = 'U (m/s)'
if field_num > 1
   start_var = start_var + bt*(we + 1)*sn
   nx = we
   ny = sn + 1
var_units = 'V (m/s)'
end
if field_num > 2
   start_var = start_var + bt*we*(sn + 1)
   nx = we
   ny = sn
var_units = 'W (m/s)'
end
if field_num > 3
   start_var = start_var + (bt + 1)*we*sn
var_units = 'GZ (m^2/s^2)'
end
if field_num > 4
   start_var = start_var + (bt + 1)*we*sn
var_units = 'T (K)'
end
if field_num > 5
   start_var = start_var + bt*we*sn
var_units = 'MU (Pa)'
end
if field_num > 6
   start_var = start_var + we*sn
var_units = 'QV (kg/kg)'
end
if field_num > 7
   start_var = start_var + bt*we*sn*(field_num-7)
var_units = 'QC (kg/kg)'
end
if field_num > 8
var_units = 'QR (kg/kg)'
end

plot_title = ['True state   ' var_units '   Level: ' num2str(field_level) '   Time: ' num2str(itime)]

start_var = start_var + nx*ny*(field_level - 1)

% Extract field

field_vec = single_state(start_var : start_var + nx*ny - 1);

field = reshape(field_vec, [nx, ny]);

% Plot field

%nc=5

%colormap = (prism(nc))
[C, h] = contourf(field');
title(plot_title)
colorbar('vert')
%clabel(C, h);

% Loop for another try
%map_wrf;

