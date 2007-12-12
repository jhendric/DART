function fit_ens_mean_time(ddir)
% fit_ens_mean_time(ddir)
%
% Part of the observation-space diagnostics routines.
%
% Plots the spatial mean RMSE of the ensemble mean as a function of 
% time for both the 'guess' and the 'analysis' at a single level.
% Several regions are plotted. This function simply plots the
% data in *ges_times.dat using metadata in ObsDiagAtts.m - both
% created by the executable 'obs_diag'.
%
% 'obs_diag' also produces a matlab-compatible file of plotting attributes:
% ObsDiagAtts.m which specifies the run-time configuration of obs_diag.
%
% ddir     is an optional argument specifying the directory containing
%               the data files as preprocessed by the support routines.
%
% USAGE: if the preprocessed data files are in a directory called 'plot'
%
% ddir = 'plot';
% fit_ens_mean_time(ddir)
%
% USAGE: if the preprocessed data files are in the current directory 
%
% fit_ens_mean_time

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% Ensures the specified directory is searched first.
if ( nargin > 0 )
   startpath = addpath(ddir);
else
   startpath = path;
end

%----------------------------------------------------------------------
% Defaults
%----------------------------------------------------------------------

datafile = 'ObsDiagAtts';
ptypes = {'gs-','bd-','ro-','k+-'};    % for each region

%----------------------------------------------------------------------
% Get plotting metadata from obs_diag run.
%----------------------------------------------------------------------

if ( exist(datafile) == 2 )

   eval(datafile)

   if ( exist('plevel','var') == 0 )
      plevel = 1;
      iskip = iskip_days;
      plotdat.toff = 0;
      plotdat.bin1 = datenum(t1);
   else  % high dimensional models
      % Coordinate between time types and dates
      skip_seconds = time_to_skip(4)*3600 + time_to_skip(5)*60 + time_to_skip(6);
      iskip = time_to_skip(3) + skip_seconds/86400;
      plotdat.bin1 = datenum(first_bin_center); % a known date in matlab's time units
      plotdat.toff = plotdat.bin1 - t1;
   end

else
   error(sprintf('%s cannot be found.', datafile))
end

% Set up a structure with all the plotting components
plotdat.day1      = datestr(t1+plotdat.toff+iskip,'yyyy-mm-dd HH');
plotdat.dayN      = datestr(tN+plotdat.toff,'yyyy-mm-dd HH');
plotdat.level     = plevel;
plotdat.ylabel    = 'RMSE';
plotdat.nregions  = length(Regions);
plotdat.nvars     = length(One_Level_Varnames);
plotdat.flavor    = 'Ens Mean';

%----------------------------------------------------------------------
% Loop around observation types
%----------------------------------------------------------------------

for ivar = 1:plotdat.nvars,

   plotdat.varname = One_Level_Varnames{ivar};

   switch obs_select
      case 1,
         string1 = sprintf('%s (all data)',     plotdat.varname);
      case 2, 
         string1 = sprintf('%s (RaObs)',        plotdat.varname);
      otherwise,
         string1 = sprintf('%s (ACARS,SATWND)', plotdat.varname);
   end

%  switch One_Level_Varnames{ivar}
%     case{'P'}
         ges  = sprintf('%s_ges_times.dat',One_Level_Varnames{ivar});
         anl  = sprintf('%s_anl_times.dat',One_Level_Varnames{ivar});
         main = sprintf('%s %s',plotdat.flavor,string1);
%     otherwise
%        ges  = sprintf('%s_ges_times_%04dmb.dat',One_Level_Varnames{ivar});
%        anl  = sprintf('%s_anl_times_%04dmb.dat',One_Level_Varnames{ivar});
%        main = sprintf('%s %s %d hPa',plotdat.flavor,string1,plotdat.level);
%  end

   plotdat.ges     = ges;
   plotdat.anl     = anl;

   % plot each region

   figure(ivar); clf; 

   for iregion = 1:length(Regions),
      plotdat.title  = Regions{iregion};
      plotdat.region = iregion;
      myplot(plotdat);
   end

   CenterAnnotation(main);  % One title in the middle
   BottomAnnotation(ges);   % annotate filename at bottom

   % create a postscript file

   psfname = sprintf('%s_ens_mean_time.ps',plotdat.varname);
   print(ivar,'-dpsc',psfname);

end

path(startpath); % restore MATLABPATH to original setting

%----------------------------------------------------------------------
% 'Helper' functions
%----------------------------------------------------------------------

function myplot(plotdat)
p1 = load(plotdat.ges); p = SqueezeMissing(p1); 
a1 = load(plotdat.anl); a = SqueezeMissing(a1); 

% Make x axis plotting arrays in units of 'days'
xp = p(:,1) + p(:,2)/86400 + plotdat.toff;
xa = a(:,1) + a(:,2)/86400 + plotdat.toff;

offset = 3;  % columns 1,2 are time, 3=mean, 4=spread, 5=numobs

count  = offset+(plotdat.region-1)*3;
yp     = p(:,count);
ya     = a(:,count);

gmean = mean(yp(isfinite(yp))); gstring = sprintf('guess;    mean=%.3f',gmean);
amean = mean(ya(isfinite(ya))); astring = sprintf('analysis; mean=%.3f',amean);

if ( plotdat.nregions > 2 )
   subplot(2,2,plotdat.region)
else
   subplot(plotdat.nregions,1,plotdat.region)
end

   plot(xp,yp,'k+-',xa,ya,'ro-','LineWidth',1.5)
   grid
   ax = axis; ax(3) = 0.0; axis(ax);
   ylabel(plotdat.ylabel, 'fontsize', 10);
   title(plotdat.title, 'Interpreter','none','fontsize', 12,'FontWeight','bold')
   h = legend(gstring, astring);
   legend(h,'boxoff')

   % a slightly better way to annotate dates, etc.
   ttot = max(xp) - min(xp) + 1;
   if ((plotdat.bin1 > 1000) && (ttot > 32));
      datetick('x',6,'keeplimits','keepticks');
      monstr = datestr(xp(1),28);
      xlabel(sprintf('month/day - %s start',monstr))
   elseif (plotdat.bin1 > 1000);
      datetick('x',7);
      monstr = datestr(xp(1),28);
      xlabel(sprintf('day of month - %s start',monstr))
   else
      xlabel('days')
   end



function y = SqueezeMissing(x)

missing = find(x < -98); % 'missing' is coded as -99

if isempty(missing)
  y = x;
else
  y = x;
  y(missing) = NaN;
end



function CenterAnnotation(main)
subplot('position',[0.48 0.48 0.04 0.04])
axis off
h = text(0.5,0.5,main);
set(h,'HorizontalAlignment','center', ...
      'VerticalAlignment','bottom', ...
      'Interpreter','none', ...
      'FontSize',12, ...
      'FontWeight','bold')



function BottomAnnotation(main)
% annotates the directory containing the data being plotted
subplot('position',[0.48 0.01 0.04 0.04])
axis off
bob = which(main);
[pathstr,name,ext,versn] = fileparts(bob);
h = text(0.0,0.5,pathstr);
set(h,'HorizontalAlignment','center', ...
      'VerticalAlignment','middle',...
      'Interpreter','none',...
      'FontSize',8)