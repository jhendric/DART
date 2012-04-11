function PlotCorrel( pinfo )
%% PlotCorrel   space-time series of correlation between a variable at a given
% time and all variables at all times in an ensemble time sequence.
%
% PlotCorrel is intended to be called by 'plot_correl'.
%
% USAGE: PlotCorrel( pinfo )
%
% pinfo      A structure containing all necessary plotting information.
%            For the low-order models, the structure MUST contain:
%
% fname             name of netCDF file containing a DART ensemble
% base_var          name of netCDF variable
% base_var_index    index of state variable used as standard in correlation
% base_time         index of time series to use as the standard for correlation
%
% Example 1   (9var model with 1000 time steps)
%%------------------------------------------------------------------
% pinfo.fname          = 'Prior_Diag.nc';
% pinfo.base_var       = 'state';
% pinfo.base_var_index = 5;          % picked arbitrarily
% pinfo.base_time      = 238;        % ditto
% PlotCorrel(pinfo)                  % generates a plot

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if (exist(pinfo.fname,'file') ~= 2), error('%s does not exist.',pinfo.fname), end

contourlevels = [-1:0.2:-0.2 0.2:0.2:1.0];  % no contour at zero, please

switch(lower(pinfo.model))

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
	 'lorenz_04','forced_lorenz_96','ikeda','simple_advection'}

      % The Base Variable Index must be a valid state variable
      if ( pinfo.base_var_index > pinfo.num_state_vars )
         fprintf('%s only has %d state variables\n', pinfo.fname, pinfo.num_state_vars)
         error('you wanted variable # %d ', pinfo.base_var_index)
      end

      % The Time must be within range also.
      if ( pinfo.base_time > pinfo.time_series_length )
         fprintf('%s only has %d output times\n', pinfo.fname, pinfo.time_series_length)
         error('you wanted time # %d ', pinfo.base_time)
      end

      % Get 'standard' ensemble series
      base = get_hyperslab('fname',pinfo.fname, 'varname',pinfo.base_var, ...
                  'copyindex1',pinfo.ensemble_indices(1),'copycount',pinfo.num_ens_members, ...
                  'stateindex',pinfo.base_var_index);

      % Get (potentially large) state.
      state_var = get_hyperslab('fname',pinfo.fname, 'varname',pinfo.base_var, ...
                  'copyindex1',pinfo.ensemble_indices(1),'copycount',pinfo.num_ens_members, ...
                  'state1',pinfo.min_state_var,'statecount',pinfo.num_state_vars);
      
      % It is efficient to preallocate correl storage ...
      correl = zeros(pinfo.num_state_vars,pinfo.time_series_length);

      % Need to loop through all variables in the ensemble
      for i = 1:pinfo.num_state_vars,
         correl(i, :) = ens_correl(base, pinfo.base_time, state_var(:,:,i));
      end

      % Now for the plotting part ...
      disp('Please be patient ... this usually takes a bit ...')
      clf;

      [cs,h] = contour(pinfo.time, 1:pinfo.num_state_vars, correl, contourlevels);
 %    clabel(cs,h,'FontSize',12,'Color','k','Rotation',0);
      set(gca,'Clim',[-1 1])
      hold on; % highlight the reference state variable and time
      plot(pinfo.time(pinfo.base_time), pinfo.base_var_index,'kh','MarkerSize',12,'MarkerFaceColor','k')

      s1 = sprintf('%s Correlation of variable %s index %d, timestep = %d', ...
               pinfo.model, pinfo.base_var, pinfo.base_var_index, pinfo.base_time);
      s2 = sprintf('against all variables, all times, %d ensemble members', ...
               size(state_var,2));
      title({s1,s2,pinfo.fname},'interpreter','none','fontweight','bold')
      xlabel(sprintf('model "days" (%d timesteps)',pinfo.time_series_length))
      ylabel('state variable (index)')
      set(gca,'YTick',1:pinfo.num_state_vars)
      colorbar


   case 'fms_bgrid'

      % We are going to correlate one var/time/lvl/lat/lon  with
      % all other lats/lons for a var/time/lvl

      clf;

      switch lower(pinfo.comp_var)
         case {'ps','t'}
            lats     = nc_varget(pinfo.fname,'TmpJ'); ny = length(lats);
            lons     = nc_varget(pinfo.fname,'TmpI'); nx = length(lons);
            latunits = nc_attget(pinfo.fname,'TmpJ','units');
            lonunits = nc_attget(pinfo.fname,'TmpI','units');
         otherwise
            lats     = nc_varget(pinfo.fname,'VelJ'); ny = length(lats);
            lons     = nc_varget(pinfo.fname,'VelI'); nx = length(lons);
            latunits = nc_attget(pinfo.fname,'VelJ','units');
            lonunits = nc_attget(pinfo.fname,'VelI','units');
      end

      nxny = nx*ny;

      base_mem = get_hyperslab('fname', pinfo.fname, 'varname', pinfo.base_var, ...
                    'copyindex1', pinfo.ensemble_indices(1), 'copycount',pinfo.num_ens_members, ...
                    'timeindex', pinfo.base_tmeind, 'levelindex',pinfo.base_lvlind, ...
                    'latindex',  pinfo.base_latind, 'lonindex',  pinfo.base_lonind );

      bob      = get_hyperslab('fname', pinfo.fname, 'varname', pinfo.comp_var, ...
                    'copyindex1', pinfo.ensemble_indices(1), 'copycount',pinfo.num_ens_members, ...
                    'timeindex', pinfo.base_tmeind, 'levelindex', pinfo.comp_lvlind);

      comp_ens = reshape(bob,[pinfo.num_ens_members, nxny]);
      corr     = zeros(nxny,1);

      for i = 1:nxny,
         x = corrcoef(base_mem, comp_ens(:, i));
         corr(i) = x(1, 2);
      end

      correl = reshape(corr,[ny nx]);

      [cs, h] = contour(lons,lats,correl,contourlevels);
      clabel(cs,h,'FontSize',12,'Color','k','Rotation',0);
      set(gca,'Clim',[-1 1])
      hold on;
      plot(pinfo.base_lon, pinfo.base_lat, 'pk', ...
                 'MarkerSize',12,'MarkerFaceColor','k');
      s1 = sprintf('%s Correlation of ''%s'', level %d, (%.2f,%.2f) T = %f', ...
           pinfo.model, pinfo.base_var, pinfo.base_lvl, ...
             pinfo.base_lat, pinfo.base_lon, pinfo.base_time);

      s2 = sprintf('against ''%s'', entire level %d, same time, %d ensemble members', ...
               pinfo.comp_var, pinfo.comp_lvl, pinfo.num_ens_members);
      title({s1,s2,pinfo.fname},'interpreter','none','fontweight','bold')
      xlabel(sprintf('longitude (%s)',lonunits),'interpreter','none')
      ylabel(sprintf('latitude (%s)',latunits),'interpreter','none')
      continents;
      axis image
      colorbar;

   case 'wrf'

      % We are going to correlate one var/time/lvl/lat/lon  with
      % all other lats/lons for a var/time/lvl

      clf;

      % Get the plotting lat/lon for the comparison variable.
      % This is the variable that has a spatial extent.

      varinfo = nc_getvarinfo(pinfo.fname, pinfo.comp_var);
      latdim  = find(strncmp('south_north',varinfo.Dimension,length('south_north')) > 0);
      londim  = find(strncmp(  'west_east',varinfo.Dimension,length(  'west_east')) > 0);
      ny      = varinfo.Size(latdim);
      nx      = varinfo.Size(londim);
      nxny    = nx*ny;

      % Each of the WRF variables has a 'coordinate' attribute signifying which
      % of the 6 possible lat/lon variables is appropriate.

      coordinates{1} = sscanf(nc_attget(pinfo.fname,pinfo.comp_var,'coordinates'),'%s %*s');
      coordinates{2} = sscanf(nc_attget(pinfo.fname,pinfo.comp_var,'coordinates'),'%*s %s');
      latcoord = find(strncmp('XLAT',coordinates,length('XLAT')) > 0);
      loncoord = find(strncmp('XLON',coordinates,length('XLON')) > 0);
      latmat   = nc_varget(pinfo.fname,coordinates{latcoord});
      lonmat   = nc_varget(pinfo.fname,coordinates{loncoord});
      latunits = nc_attget(pinfo.fname,coordinates{latcoord},'units');
      lonunits = nc_attget(pinfo.fname,coordinates{latcoord},'units');

      inds = (lonmat < 0); % Convert to 0,360 to minimize dateline probs.
      lonmat(inds) = lonmat(inds) + 360.0;
      if (pinfo.base_lon < 0), pinfo.base_lon = pinfo.base_lon + 360.0; end

      % Get the actual goods ... and perform the correlation

      base_mem = get_hyperslab('fname', pinfo.fname, 'varname', pinfo.base_var, ...
                    'copyindex1', pinfo.ensemble_indices(1), 'copycount',pinfo.num_ens_members, ...
                    'timeindex', pinfo.base_tmeind, 'levelindex',pinfo.base_lvlind, ...
                    'latindex',  pinfo.base_latind, 'lonindex',  pinfo.base_lonind );

      if (std(base_mem) == 0.0)
          warning('%s at level %d lat %d lon %d time %s is a constant\n',pinfo.base_var,...
             pinfo.base_lvlind,pinfo.base_latind,pinfo.base_lonind,datestr(pinfo.base_time))
          error('Cannot calculate correlation coefficient with a constant.')
      end

      bob      = get_hyperslab('fname', pinfo.fname, 'varname', pinfo.comp_var, ...
                    'copyindex1', pinfo.ensemble_indices(1), 'copycount',pinfo.num_ens_members, ...
                    'timeindex', pinfo.base_tmeind, 'levelindex', pinfo.comp_lvlind);

      if (std(bob(:)) == 0.0)
          warning('%s at level %d time %s is a constant\n',pinfo.comp_var,...
             pinfo.comp_lvlind, datestr(pinfo.base_time))
          error('Cannot calculate correlation coefficient with a constant.')
      end

      comp_ens = reshape(bob,[pinfo.num_ens_members, nxny]);
      corr     = zeros(nxny,1);

      % Really should check to see if each comp_ens is a constant value as
      % well - this is slow enough already.

      fprintf('Performing correlations at %d locations ...\n',nxny)
      for i = 1:nxny,
         x = corrcoef(base_mem, comp_ens(:, i));
         corr(i) = x(1, 2);
      end

      correl = reshape(corr,[ny nx]);

      % Plot it up ...

      [cs,h] = contour(lonmat,latmat,correl,contourlevels);
      clabel(cs,h,'FontSize',12,'Color','k','Rotation',0);
      set(gca,'Clim',[-1 1])
      hold on;
      plot(pinfo.base_lon, pinfo.base_lat, 'pk', ...
                 'MarkerSize',12,'MarkerFaceColor','k');
      s1 = sprintf('%s Correlation of ''%s'', level %d, (%.2f,%.2f) T = %s', ...
           pinfo.model, pinfo.base_var, pinfo.base_lvl, ...
             pinfo.base_lat, pinfo.base_lon, datestr(pinfo.base_time));

      s2 = sprintf('against ''%s'', entire level %d, same time, %d ensemble members', ...
               pinfo.comp_var, pinfo.comp_lvl, pinfo.num_ens_members);
      title({s1,s2,pinfo.fname},'interpreter','none','fontweight','bold')
      xlabel(sprintf('longitude (%s)',lonunits),'interpreter','none')
      ylabel(sprintf('latitude (%s)',latunits),'interpreter','none')
      continents;
      axis image
      colorbar;

   case 'mitgcm_ocean'

      % We are going to correlate one var/time/lvl/lat/lon  with
      % all other lats/lons for a var/time/lvl

      clf;

      switch lower(pinfo.comp_var)
         case {'u'}
            lats     = nc_varget(pinfo.fname,'YC'); ny = length(lats);
            lons     = nc_varget(pinfo.fname,'XG'); nx = length(lons);
            latunits = nc_attget(pinfo.fname,'YC','units');
            lonunits = nc_attget(pinfo.fname,'XG','units');
         case {'v'}
            lats     = nc_varget(pinfo.fname,'YG'); ny = length(lats);
            lons     = nc_varget(pinfo.fname,'XC'); nx = length(lons);
            latunits = nc_attget(pinfo.fname,'YG','units');
            lonunits = nc_attget(pinfo.fname,'XC','units');
         otherwise
            lats     = nc_varget(pinfo.fname,'YC'); ny = length(lats);
            lons     = nc_varget(pinfo.fname,'XC'); nx = length(lons);
            latunits = nc_attget(pinfo.fname,'YC','units');
            lonunits = nc_attget(pinfo.fname,'XC','units');
      end

      nxny = nx*ny;

      base_mem = get_hyperslab('fname', pinfo.fname, 'varname', pinfo.base_var, ...
                    'copyindex1', pinfo.ensemble_indices(1), 'copycount',pinfo.num_ens_members, ...
                    'timeindex', pinfo.base_tmeind, 'levelindex',pinfo.base_lvlind, ...
                    'latindex',  pinfo.base_latind, 'lonindex',  pinfo.base_lonind );

      bob      = get_hyperslab('fname', pinfo.fname, 'varname', pinfo.comp_var, ...
                    'copyindex1', pinfo.ensemble_indices(1), 'copycount',pinfo.num_ens_members, ...
                    'timeindex', pinfo.base_tmeind, 'levelindex', pinfo.comp_lvlind);

      comp_ens = reshape(bob,[pinfo.num_ens_members, nxny]);
      corr     = zeros(nxny,1);

      for i = 1:nxny,
         x = corrcoef(base_mem, comp_ens(:, i));
         corr(i) = x(1, 2);
      end

      correl = reshape(corr,[ny nx]);

      [cs,h] = contour(lons,lats,correl,contourlevels);
      clabel(cs,h,'FontSize',12,'Color','k','Rotation',0);
      set(gca,'Clim',[-1 1])
      hold on;
      plot(pinfo.base_lon, pinfo.base_lat, 'pk', ...
                 'MarkerSize',12,'MarkerFaceColor','k');
      s1 = sprintf('%s Correlation of ''%s'', level %d, (%.2f,%.2f) T = %f', ...
           pinfo.model, pinfo.base_var, pinfo.base_lvl, ...
             pinfo.base_lat, pinfo.base_lon, pinfo.base_time);

      s2 = sprintf('against ''%s'', entire level %d, same time, %d ensemble members', ...
               pinfo.comp_var, pinfo.comp_lvl, pinfo.num_ens_members);
      title({s1,s2,pinfo.fname},'interpreter','none','fontweight','bold')
      xlabel(sprintf('longitude (%s)',lonunits),'interpreter','none')
      ylabel(sprintf('latitude (%s)',latunits),'interpreter','none')
      continents;
      axis image
      colorbar;

   case {'pe2lyr','cam'}

      % We are going to correlate one var/time/lvl/lat/lon  with
      % all other lats/lons for a var/time/lvl

      clf;

      lats     = nc_varget(pinfo.fname,'lat'); ny = length(lats);
      lons     = nc_varget(pinfo.fname,'lon'); nx = length(lons);
      latunits = nc_attget(pinfo.fname,'lat','units');
      lonunits = nc_attget(pinfo.fname,'lon','units');
      nxny     = nx*ny;

      base_mem = get_hyperslab('fname', pinfo.fname, 'varname', pinfo.base_var, ...
                    'copyindex1', pinfo.ensemble_indices(1), 'copycount',pinfo.num_ens_members, ...
                    'timeindex', pinfo.base_tmeind, 'levelindex',pinfo.base_lvlind, ...
                    'latindex',  pinfo.base_latind, 'lonindex',  pinfo.base_lonind );

      bob      = get_hyperslab( 'fname', pinfo.fname, 'varname', pinfo.comp_var, ...
                    'copyindex1', pinfo.ensemble_indices(1), 'copycount',pinfo.num_ens_members, ...
                    'timeindex', pinfo.base_tmeind, 'levelindex', pinfo.comp_lvlind);

      comp_ens = reshape(bob,[pinfo.num_ens_members, nxny]);
      corr     = zeros(nxny,1);

      for i = 1:nxny,
         x = corrcoef(base_mem, comp_ens(:, i));
         corr(i) = x(1, 2);
      end

      correl = reshape(corr,[ny nx]);

      [cs,h] = contour(lons,lats,correl,contourlevels);
      clabel(cs,h,'FontSize',12,'Color','k','Rotation',0);
      set(gca,'Clim',[-1 1])
      hold on;
      plot(pinfo.base_lon, pinfo.base_lat, 'pk', ...
                 'MarkerSize',12,'MarkerFaceColor','k');
      s1 = sprintf('%s Correlation of ''%s'', level %d, (%.2f,%.2f) T = %f', ...
           pinfo.model, pinfo.base_var, pinfo.base_lvl, ...
             pinfo.base_lat, pinfo.base_lon, pinfo.base_time);

      s2 = sprintf('against ''%s'', entire level %d, same time, %d ensemble members', ...
               pinfo.comp_var, pinfo.comp_lvl, pinfo.num_ens_members);
      title({s1,s2,pinfo.fname},'interpreter','none','fontweight','bold')
      xlabel(sprintf('longitude (%s)',lonunits),'interpreter','none')
      ylabel(sprintf('latitude (%s)',latunits),'interpreter','none')
      continents;
      axis image
      colorbar;

   case {'tiegcm'}

      % We are going to correlate one var/time/lvl/lat/lon  with
      % all other lats/lons for a var/time/lvl

      clf;

      lats     = nc_varget(pinfo.fname,'lat'); ny = length(lats);
      lons     = nc_varget(pinfo.fname,'lon'); nx = length(lons);
      latunits = nc_attget(pinfo.fname,'lat','units');
      lonunits = nc_attget(pinfo.fname,'lon','units');

      inds = find(lons >= 180);
      lons(inds) = lons(inds) - 360.0;

      nxny     = nx*ny;

      base_mem = get_hyperslab('fname', pinfo.fname, 'varname', pinfo.base_var, ...
                    'copyindex1', pinfo.ensemble_indices(1), 'copycount',pinfo.num_ens_members, ...
                    'timeindex', pinfo.base_tmeind, 'levelindex',pinfo.base_lvlind, ...
                    'latindex',  pinfo.base_latind, 'lonindex',  pinfo.base_lonind );

      bob      = get_hyperslab('fname', pinfo.fname, 'varname', pinfo.comp_var, ...
                    'copyindex1', pinfo.ensemble_indices(1), 'copycount',pinfo.num_ens_members, ...
                    'timeindex', pinfo.base_tmeind, 'levelindex', pinfo.comp_lvlind);

      comp_ens = reshape(bob,[pinfo.num_ens_members, nxny]);
      corr     = zeros(nxny,1);

      for i = 1:nxny,
         x = corrcoef(base_mem, comp_ens(:, i));
         corr(i) = x(1, 2);
      end

      correl = reshape(corr,[ny nx]);

%     imagesc(lons,lats,correl); set(gca,'YDir','normal'); hold on;
      [cs,h] = contour(lons,lats,correl,contourlevels);
      clabel(cs,h,'FontSize',12,'Color','k','Rotation',0);
      set(gca,'Clim',[-1 1])
      hold on;
      plot(pinfo.base_lon, pinfo.base_lat, 'pk', ...
                 'MarkerSize',12,'MarkerFaceColor','k');
      s1 = sprintf('%s Correlation of ''%s'', level %d, (%.2f,%.2f) T = %f', ...
           pinfo.model, pinfo.base_var, pinfo.base_lvl, ...
             pinfo.base_lat, pinfo.base_lon, pinfo.base_time);

      s2 = sprintf('against ''%s'', entire level %d, same time, %d ensemble members', ...
               pinfo.comp_var, pinfo.comp_lvl, pinfo.num_ens_members);
      title({s1,s2,pinfo.fname},'interpreter','none','fontweight','bold')
      xlabel(sprintf('longitude (%s)',lonunits),'interpreter','none')
      ylabel(sprintf('latitude (%s)',latunits),'interpreter','none')
      continents('hollow','dateline');
      axis image
      colorbar;

   case {'mpas_atm'}

      %% We are going to correlate one var/time/lvl/location  with
      %  all other locations for a var/time/lvl

      clf;

      base_mem = get_hyperslab('fname',pinfo.fname, 'varname',pinfo.base_var, ...
                    'copyindex1', pinfo.ensemble_indices(1), 'copycount',pinfo.num_ens_members, ...
                    'timeindex', pinfo.base_tmeind, 'levelindex',pinfo.base_lvlind, ...
                    'cellindex', pinfo.base_cellindex );

      comp_ens = get_hyperslab('fname',pinfo.fname, 'varname',pinfo.comp_var, ...
                    'copyindex1', pinfo.ensemble_indices(1), 'copycount',pinfo.num_ens_members, ...
                    'timeindex', pinfo.base_tmeind, 'levelindex', pinfo.comp_lvlind);

      [~, nxny] = size(comp_ens);
      corr      = zeros(nxny,1);

      for i = 1:nxny,
         x = corrcoef(base_mem, comp_ens(:, i));
         corr(i) = x(1, 2);
      end

      %% here's the tricky part ... plotting the unstructured grid

      PlotMPAScells(pinfo.fname, corr)
      hold on
      plot(pinfo.lonCell(pinfo.base_cellindex), pinfo.latCell(pinfo.base_cellindex),'pb','MarkerSize',20);
      hold off

      s1 = sprintf('%s Correlation of ''%s'', level %d, (%.2f,%.2f) T = %s', ...
           pinfo.model, pinfo.base_var, pinfo.base_lvl, ...
           pinfo.latCell(pinfo.base_cellindex), pinfo.lonCell(pinfo.base_cellindex), ...
           datestr(pinfo.base_time));

      s2 = sprintf('against ''%s'', entire level %d, same time, %d ensemble members', ...
               pinfo.comp_var, pinfo.comp_lvl, pinfo.num_ens_members);
      title({s1,s2,pinfo.fname},'interpreter','none','fontweight','bold')
      xlabel(sprintf('longitude (%s)',pinfo.lonunits),'interpreter','none')
      ylabel(sprintf('latitude (%s)',pinfo.latunits),'interpreter','none')
      continents('hollow');
      axis image
      set(gca,'Clim',[-1 1])
      axis([-10 370 -Inf Inf])
      colorbar;

   otherwise

      error('model %s not implemented yet', pinfo.model)

end


%----------------------------------------------------------------------
% helper functions
%----------------------------------------------------------------------


function PlotMPAScells(fname,x)

%% Get the plotting arrays
%  for some reason, the lonCell,latCell stuff is [0,360] and the
%  latVertex,lonVertex is [-pi,pi] ... daffy.

nEdgesOnCell   = nc_varget(fname,'nEdgesOnCell');
verticesOnCell = nc_varget(fname,'verticesOnCell');
latVertex      = nc_varget(fname,'latVertex') * 180/pi; % cvrt to [-90,90]
lonVertex      = nc_varget(fname,'lonVertex') * 180/pi; % cvrt to [-180,180]

inds = lonVertex < 0;
lonVertex(inds) = lonVertex(inds) + 360.0;

%% Each Cell has some number of vertices

nrows  = max( nEdgesOnCell);
nCells = size(nEdgesOnCell,1);
xpoly  = NaN(nrows,nCells);
ypoly  = NaN(nrows,nCells);

for iCell=1:nCells

   n = nEdgesOnCell(iCell);

   for i = 1:n
      xpoly(i,iCell) = lonVertex(verticesOnCell(iCell,i));
      ypoly(i,iCell) = latVertex(verticesOnCell(iCell,i));
      if (i > 1)
         if (abs(xpoly(i,iCell) - xpoly(1,iCell)) > 180.0)
            if (xpoly(i,iCell) > xpoly(1,iCell))
                xpoly(i,iCell) = xpoly(i,iCell) - 360.0;
            else
                xpoly(i,iCell) = xpoly(i,iCell) + 360.0;
            end
         end
      end
   end

   % patch the pentagons up to hexagons
   if (n < nrows)
      xpoly(n+1:nrows,iCell) = xpoly(n,iCell);
      ypoly(n+1:nrows,iCell) = ypoly(n,iCell);
   end

end

h = patch(xpoly,ypoly,x');
set(h,'LineStyle','none');

