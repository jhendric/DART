% psfc_movie
%

% Data Assimilation Research Testbed -- DART
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

field_name = 'MU';

fname = 'psfc';

nc = netcdf( [fname,'.nc'] , 'nowrite' ) ;

we = size(nc('west_east'),1);
sn = size(nc('south_north'),1);
dt = nc.DT(:);

close(nc);

f_size = we*sn;

true_times = getnc(fname, 'time');
num_true_times = size(true_times, 1)

stime = input('Initial time : ');
ftime = input('End time : ');

var_units = ' (Pa)';
iso = [-2:0.1:2];

%set(gca,'nextplot','replacechildren');

pane = 1;

for itime = stime:ftime

% Extract field

   field1 = getnc(fname, field_name,[itime -1 -1],[itime -1 -1],[1 1 1]);
   field2 = getnc(fname, field_name,[itime+1 -1 -1],[itime+1 -1 -1],[1 1 1]);

   field = (field2 - field1)/dt;

% Plot field

   if min(min(field)) ~= max(max(field))

      [C h] = contourf(field, iso); caxis([min(iso(:)),max(iso(:))]);
%     h = clabel(C,h,'labelspacing',288);  set(h,'Fontsize',12);
%     hold on
%     [Cm hm] = contourm(xlat,xlon,field, -iso, 'b--','LineWidth',2);
%     hm = clabelm(Cm,hm,'labelspacing',288);  set(hm,'Fontsize',12);

   plot_title = [ field_name var_units ...
			    '   ' true_times(itime,:) ];

   title(plot_title)
   colorbar

   F(pane) = getframe(gcf);

   pane = pane + 1;

   end

end
