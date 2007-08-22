function PlotTotalErr( pinfo )
% PlotTotalErr Plots summary plots of global error and spread
%
% The error is a measure of the distance between the 
% 'ensemble mean' and the 'true_state'. The distance of the 
% spread is also plotted.
%
% PlotTotalErr is intended to be called by 'plot_total_err'
% The only input argument is a structure with model-dependent
% components.
%
% USAGE: PlotTotalErr( pinfo )
%
% STRUCTURE COMPONENTS FOR low-order models
% truth_file      name of netCDF DART file with copy tagged 'true state'
% diagn_file      name of netCDF DART file with copies tagged 'ensemble mean'
%                 and 'ensemble spread'
%
% Example 1   (Lorenz_63, Lorenz_96, 9var ...)
%%--------------------------------------------------------
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Posterior_Diag.nc';
% PlotTotalErr( pinfo )

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

% adds the time overlap info to the structure.  
pstruct = CheckModelCompatibility(pinfo);
pinfo = CombineStructs(pinfo, pstruct) 

f = netcdf(pinfo.truth_file,'nowrite');
model = f.model(:); 

% Get the netcdf variable indices for desired "copies"
% Get the indices for the true state, ensemble mean and spread
% The metadata is queried to determine which "copy" is appropriate.
truth_index      = get_copy_index(pinfo.truth_file, 'true state' );
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% The times are set before from the model compatibility call
times = getnc(pinfo.truth_file,'time', [pinfo.truth_time(1)], ...
                                       [pinfo.truth_time(2)]);
num_times = length(times);

switch lower(model)

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_04','ikeda'}

      % Get the appropriate netcdf variables
      truth  = get_state_copy(pinfo.truth_file, 'state',     truth_index, ...
                              pinfo.truth_time(1), pinfo.truth_time(2));
      ens    = get_state_copy(pinfo.diagn_file, 'state',  ens_mean_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
      spread = get_state_copy(pinfo.diagn_file, 'state',ens_spread_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
      num_vars = size(spread,2);

      % Also need to compute the spread; zero truth for this and
      % compute distance from 0
      err        = total_err(truth, ens);
      err_spread = total_err(zeros(size(spread)), spread);
      errTotal   = sum(err)/num_times;
      spreadTotal= sum(err_spread)/num_times;
      string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
      string2 = ['time-mean Ensemble Spread Total Error = ' num2str(spreadTotal)];

      clf;
      plot(times,err, 'b', times,err_spread, 'r');
      legend(string1,string2,0)
      legend boxoff
      title(sprintf('%s Total Error over all %d variables for %s',...
                    model, num_vars, pinfo.diagn_file), ...
            'interpreter','none','fontweight','bold')
      xlabel(sprintf('model time (%d timesteps)',num_times))
      ylabel('Total Error')

   case 'lorenz_96_2scale'

      % Simply going to append X,Y together and treat as above.

      % Get the appropriate netcdf variables
      tim    = get_state_copy(pinfo.truth_file, 'X',     truth_index, ...
                              pinfo.truth_time(1), pinfo.truth_time(2));
      tom    = get_state_copy(pinfo.truth_file, 'Y',     truth_index, ...
                              pinfo.truth_time(1), pinfo.truth_time(2));
      truth  = [tim tom];
      tim    = get_state_copy(pinfo.diagn_file, 'X',  ens_mean_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
      tom    = get_state_copy(pinfo.diagn_file, 'Y',  ens_mean_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
      ens    = [tim tom]; 
      tim    = get_state_copy(pinfo.diagn_file, 'X',ens_spread_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
      tom    = get_state_copy(pinfo.diagn_file, 'Y',ens_spread_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
      spread = [tim tom]; clear tim tom
      num_vars = size(spread,2);

      % Also need to compute the spread; zero truth for this and
      % compute distance from 0
      err        = total_err(truth, ens);
      err_spread = total_err(zeros(size(spread)), spread);
      errTotal   = sum(err)/num_times;
      spreadTotal= sum(err_spread)/num_times;
      string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
      string2 = ['time-mean Ensemble Spread Total Error = ' num2str(spreadTotal)];

      clf;
      plot(times,err, 'b', times,err_spread, 'r');
      legend(string1,string2,0)
      legend boxoff
      title(sprintf('%s Total Error over all %d variables for %s',...
                    model, num_vars, pinfo.diagn_file), ...
            'interpreter','none','fontweight','bold')
      xlabel(sprintf('model time (%d timesteps)',num_times))
      ylabel('Total Error')

   case 'forced_lorenz_96'

      % This model has the state variables replicated, so there is a difference
      % between num_state_vars and the length of the state variable.
      forcing           = f.model_forcing(:);
      delta_t           = f.model_delta_t(:);
      time_step_days    = f.model_time_step_days(:);
      time_step_seconds = f.model_time_step_seconds(:);
      num_model_vars    = f.model_num_state_vars(:);  % ACTUAL state vars, not

      % Get the appropriate netcdf variables

      Whole_truth  = get_state_copy(pinfo.truth_file, 'state',     truth_index, ...
                              pinfo.truth_time(1), pinfo.truth_time(2));
      Whole_ens    = get_state_copy(pinfo.diagn_file, 'state',  ens_mean_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
      Whole_spread = get_state_copy(pinfo.diagn_file, 'state',ens_spread_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
      num_vars = size(Whole_spread,2);

      %--------------------------------------------------------------------------
      % Treat the traditional state variable independent of the forcing variables
      %--------------------------------------------------------------------------

      ind1 = 1;                 % ASSUMPTION: traditional state is first
      indN = num_model_vars;

      truth  = Whole_truth(  :, ind1:indN );
      ens    = Whole_ens(    :, ind1:indN );
      spread = Whole_spread( :, ind1:indN );

      % Compute the spread; zero truth for this and compute distance from 0

      err        = total_err(truth, ens);
      err_spread = total_err(zeros(size(spread)), spread);
      errTotal   = sum(err)/num_times;
      spreadTotal= sum(err_spread)/num_times;
      string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
      string2 = ['time-mean Ensemble Spread Total Error = ' num2str(spreadTotal)];

      clf; subplot(2,1,1);

      plot(times,err, 'b', times,err_spread, 'r');
      legend(string1,string2,0)
      legend boxoff
      title(sprintf('%s Total Error over statevars %d to %d for %s',...
                    model, ind1, indN, pinfo.diagn_file), ...
            'interpreter','none','fontweight','bold')
      xlabel(sprintf('model time (%d timesteps)',num_times))
      ylabel('Total Error')

      %--------------------------------------------------------------------------
      % Now for the forcing
      %--------------------------------------------------------------------------

      ind1 = num_model_vars + 1;
      indN = num_vars;

      truth  = Whole_truth(  :, ind1:indN );
      ens    = Whole_ens(    :, ind1:indN );
      spread = Whole_spread( :, ind1:indN );

      % Compute the spread; zero truth for this and compute distance from 0

      err        = total_err(truth, ens);
      err_spread = total_err(zeros(size(spread)), spread);
      errTotal   = sum(err)/num_times;
      spreadTotal= sum(err_spread)/num_times;
      string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
      string2 = ['time-mean Ensemble Spread Total Error = ' num2str(spreadTotal)];

      subplot(2,1,2)

      plot(times,err, 'b', times,err_spread, 'r');
      legend(string1,string2,0)
      legend boxoff
      title(sprintf('%s Total Error over statevars %d to %d for %s',...
                    model, ind1, indN, pinfo.diagn_file), ...
            'interpreter','none','fontweight','bold')
      xlabel(sprintf('model time (%d timesteps)',num_times))
      ylabel('Total Error')

   case {'simple_advection'}

      % if the 'state' variable exists ... then
      % 'concentration','source', and 'wind' do not.

      if ( isempty(f.state(:)))
         varlist = {'concentration','source','wind','mean_source','source_phase'};
      else
         varlist = {'state'};
      end

      for ivar = 1:length(varlist)

         % Get the appropriate netcdf variables
         truth  = get_state_copy(pinfo.truth_file, varlist{ivar},     truth_index, ...
                              pinfo.truth_time(1), pinfo.truth_time(2));
         ens    = get_state_copy(pinfo.diagn_file, varlist{ivar},  ens_mean_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
         spread = get_state_copy(pinfo.diagn_file, varlist{ivar},ens_spread_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
         num_vars = size(spread,2);

         % Also need to compute the spread; zero truth for this and
         % compute distance from 0
         err        = total_err(truth, ens);
         err_spread = total_err(zeros(size(spread)), spread);
         errTotal   = sum(err)/num_times;
         spreadTotal= sum(err_spread)/num_times;
         string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
         string2 = ['time-mean Ensemble Spread Total Error = ' num2str(spreadTotal)];

         figure(ivar); clf(ivar);
         plot(times,err, 'b', times,err_spread, 'r');
         legend(string1,string2,0)
         legend boxoff
         string1 = sprintf('%s Total Error over all %d variables', model, num_vars);
         string2 = sprintf('''%s'' %s', varlist{ivar}, pinfo.diagn_file);
         title({string1,string2},'interpreter','none','fontweight','bold')
         xlabel(sprintf('model time (%d timesteps)',num_times))
         ylabel('Total Error')

      end

   case 'fms_bgrid'

      BgridTotalError( pinfo )

   case 'pe2lyr'

      Pe2lyrTotalError( pinfo )

   case 'pbl_1d'

      PBL1DTotalError( pinfo )

   otherwise

      disp(sprintf('unknown model %s -- doing nothing',model))

end
close(f)

%=======================================================================
% helper functions
%=======================================================================

function PBL1DTotalError ( pinfo )

% Get some standard plotting arrays
 z_level = getnc(pinfo.truth_file, 'z_level',[pinfo.truth_time(1)], [pinfo.truth_time(2)]);
sl_level = getnc(pinfo.truth_file,'sl_level',[pinfo.truth_time(1)], [pinfo.truth_time(2)]);
times    = getnc(pinfo.truth_file,    'time',[pinfo.truth_time(1)], [pinfo.truth_time(2)]);
num_times  = length(times );

% Get the indices for the true state, ensemble mean and spread                  
% The metadata is queried to determine which "copy" is appropriate.             
truth_index      = get_copy_index(pinfo.truth_file, 'true state');
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% U variable

truth  = getnc(pinfo.truth_file,'U', ...
              [pinfo.truth_time(1)  truth_index -1], ...
              [pinfo.truth_time(2)  truth_index -1]);
ens    = getnc(pinfo.diagn_file,'U', ...
              [pinfo.diagn_time(1)  ens_mean_index -1], ...
              [pinfo.diagn_time(2)  ens_mean_index -1]);
spread = getnc(pinfo.diagn_file,'U', ...
              [pinfo.diagn_time(1)  ens_spread_index -1], ...
              [pinfo.diagn_time(2)  ens_spread_index -1]);

err        = total_err(              truth,    ens);
err_spread = total_err(zeros(size(spread)), spread);

y_error  = squeeze(mean(err,2));           % mean over all levels
y_spread = squeeze(mean(err_spread,2));    % mean over all levels

plot(times,y_error,'r-',times,y_spread,'g-')
title('PBL_1d mean error of U over time ... all levels.')
xlabel('days')
ylabel('mean (all levels) total error')
axis([-Inf Inf 0 Inf])



function BgridTotalError( pinfo )
%---------------------------------------------------------------------
% netcdf has no state vector, it has prognostic variables.
% We are going to plot the total error (over a horizontal slice) 
% for each variable and annotate an area-weighted total.
%---------------------------------------------------------------------

ft        = netcdf(pinfo.truth_file);
model     = ft.model(:); 
timeunits = ft{'time'}.units(:);
close(ft);

tstart = pinfo.truth_time(1);
tend   = pinfo.truth_time(2);
dstart = pinfo.diagn_time(1);
dend   = pinfo.diagn_time(2);

% Since the models are "compatible", get the info from either one.
tlons    = getnc(pinfo.truth_file, 'TmpI'); num_tlons  = length(tlons );
tlats    = getnc(pinfo.truth_file, 'TmpJ'); num_tlats  = length(tlats );
vlons    = getnc(pinfo.truth_file, 'VelI'); num_vlons  = length(vlons );
vlats    = getnc(pinfo.truth_file, 'VelJ'); num_vlats  = length(vlats );
levels   = getnc(pinfo.truth_file, 'lev' ); num_levels = length(levels);
times    = getnc(pinfo.truth_file, 'time', [tstart], [tend]); num_times  = length(times );
ens_mems = getnc(pinfo.diagn_file, 'copy'); ens_size   = length(ens_mems);

% Try to coordinate "time" ... a poor attempt, needs refining
ens_times     = getnc(pinfo.diagn_file, 'time', [dstart], [dend]); 
num_ens_times = length(ens_times);
if num_ens_times < num_times
   times     =     ens_times;
   num_times = num_ens_times;
end

% Initialize storage for error averaging
rms          = zeros(num_times, pinfo.num_state_vars, num_levels);
spread_final = zeros(num_times, pinfo.num_state_vars, num_levels);

% Get the indices for the true state, ensemble mean and spread                  
% The metadata is queried to determine which "copy" is appropriate.             
truth_index      = get_copy_index(pinfo.truth_file, 'true state');
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% Calculate weights for area-averaging.
twts = SphereWeights(tlons, tlats);   % Temperature Grid
vwts = SphereWeights(vlons, vlats);   % Velocity    Grid

% Can we afford to get the whole thing at once ???
%----------------------------------------------------------------------
% surface pressure has only one level.
% GetPS   returns a   num_times x num_lats*num_lons   2Darray.
%----------------------------------------------------------------------

disp('Processing surface pressure ...')

ivar   = 1;  % should do "ivar = find(xxx='ps')"
ilevel = 1;

truth      = GetPS(pinfo.truth_file,      truth_index, tstart, tend);
ens        = GetPS(pinfo.diagn_file,   ens_mean_index, dstart, dend);
spread     = GetPS(pinfo.diagn_file, ens_spread_index, dstart, dend);

err        = total_err(              truth,    ens, twts(:) );
err_spread = total_err(zeros(size(spread)), spread, twts(:) );

         rms(:,ivar,ilevel) = err;         % spatial mean
spread_final(:,ivar,ilevel) = err_spread;  % spatial mean

clear truth ens spread err err_spread  

%----------------------------------------------------------------------
% temperature ...  num_times x num_levels x num_lats x num_lons
% GetLevel returns a    num_times x num_lats*num_lons   2Darray.
%----------------------------------------------------------------------

for ilevel = 1:num_levels,     % Loop through all levels
for ivar=2:pinfo.num_state_vars,

   disp(sprintf('Processing level %d of %d ...',ilevel,num_levels))
   %-------------------------------------------------------------------
   % all vars organized    num_times x num_levels x num_lats x num_lons
   %-------------------------------------------------------------------

   truth  = GetLevel(pinfo.truth_file, pinfo.vars{ivar},      truth_index, ilevel, tstart, tend);
   ens    = GetLevel(pinfo.diagn_file, pinfo.vars{ivar},   ens_mean_index, ilevel, dstart, dend);
   spread = GetLevel(pinfo.diagn_file, pinfo.vars{ivar}, ens_spread_index, ilevel, dstart, dend);

   switch lower(pinfo.vars{ivar})
      case {'t'}
         err        = total_err(              truth,    ens, twts(:));
         err_spread = total_err(zeros(size(spread)), spread, twts(:));
      otherwise
         err        = total_err(              truth,    ens, vwts(:));
         err_spread = total_err(zeros(size(spread)), spread, vwts(:));
   end

            rms(:, ivar, ilevel) = err;
   spread_final(:, ivar, ilevel) = err_spread;

end
end % End of level loop

clear truth ens spread err err_spread

%----------------------------------------------------------------------
% Each variable in its own figure window
%----------------------------------------------------------------------
for ivar=1:pinfo.num_state_vars,

   figure(ivar); clf;
      ft       = netcdf(pinfo.truth_file);
      varunits = ft{pinfo.vars{ivar}}.units(:);
      close(ft);

      s1 = sprintf('%s %s Ensemble Mean for %s', model,pinfo.vars{ivar},pinfo.diagn_file);

      switch lower(pinfo.vars{ivar})
         case {'ps'}
            h1 = plot(times,          rms(:, ivar, 1), '-'); hold on;
            h2 = plot(times, spread_final(:, ivar, 1), '--');

            s{1} = sprintf('time-mean Ensemble Mean error  = %f', mean(         rms(:, ivar, 1)));
            s{2} = sprintf('time-mean Ensemble Spread = %f',      mean(spread_final(:, ivar, 1)));
         otherwise
            h1 = plot(times, squeeze(         rms(:, ivar, :)),'-'); hold on;
            h2 = plot(times, squeeze(spread_final(:, ivar, :)),'--');

            for i = 1:num_levels,
               s{i           } = sprintf('level %d error  %.3f', i,mean(         rms(:, ivar, i)));
               s{i+num_levels} = sprintf('level %d spread %.3f', i,mean(spread_final(:, ivar, i)));
            end
      end

      %h = legend([h1 h2],s); legend(h,'boxoff')
      h = legend(s); legend(h,'boxoff')
      grid on;
      xlabel(sprintf('time (%s) %d timesteps',timeunits,num_times))
      ylabel(sprintf('global-area-weighted distance (%s)',varunits))
      title(s1,'interpreter','none','fontweight','bold')

end




function Pe2lyrTotalError( pinfo )
%---------------------------------------------------------------------
% netcdf has no state vector, it has prognostic variables.
% We are going to plot the total error (over a horizontal slice) 
% for each variable and annotate an area-weighted total.
%---------------------------------------------------------------------

ft        = netcdf(pinfo.truth_file);
model     = ft.model(:); 
timeunits = ft{'time'}.units(:);
close(ft);

tstart = pinfo.truth_time(1);
tend   = pinfo.truth_time(2);
dstart = pinfo.diagn_time(1);
dend   = pinfo.diagn_time(2);

% Since the models are "compatible", get the info from either one.
tlons    = getnc(pinfo.truth_file, 'lon' ); num_tlons  = length(tlons );
tlats    = getnc(pinfo.truth_file, 'lat' ); num_tlats  = length(tlats );
levels   = getnc(pinfo.truth_file, 'lev' ); num_levels = length(levels);
times    = getnc(pinfo.truth_file, 'time', [tstart], [tend]); num_times  = length(times );
ens_mems = getnc(pinfo.diagn_file, 'copy'); ens_size   = length(ens_mems);

% Try to coordinate "time" ... a poor attempt, needs refining
ens_times     = getnc(pinfo.diagn_file, 'time', [dstart], [dend]); 
num_ens_times = length(ens_times);
if num_ens_times < num_times
   times     =     ens_times;
   num_times = num_ens_times;
end

% Initialize storage for error averaging
rms          = zeros(num_times, pinfo.num_state_vars, num_levels);
spread_final = zeros(num_times, pinfo.num_state_vars, num_levels);

% Get the indices for the true state, ensemble mean and spread                  
% The metadata is queried to determine which "copy" is appropriate.             
truth_index      = get_copy_index(pinfo.truth_file, 'true state');
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% Calculate weights for area-averaging.
twts = SphereWeights(tlons, tlats);

%----------------------------------------------------------------------
% Trying not to assume we can get the whole 3D array at once.
% GetLevel returns a    num_times x num_lats*num_lons   2Darray.
%----------------------------------------------------------------------

for ilevel = 1:num_levels,     % Loop through all levels
for ivar=1:pinfo.num_state_vars,

   disp(sprintf('Processing level %d of %d ...',ilevel,num_levels))
   %-------------------------------------------------------------------
   % all vars organized    num_times x num_levels x num_lats x num_lons
   %-------------------------------------------------------------------

   truth  = GetLevel(pinfo.truth_file, pinfo.vars{ivar},      truth_index, ilevel, tstart, tend);
   ens    = GetLevel(pinfo.diagn_file, pinfo.vars{ivar},   ens_mean_index, ilevel, dstart, dend);
   spread = GetLevel(pinfo.diagn_file, pinfo.vars{ivar}, ens_spread_index, ilevel, dstart, dend);

   err        = total_err(              truth,    ens, twts(:));
   err_spread = total_err(zeros(size(spread)), spread, twts(:));

            rms(:, ivar, ilevel) = err;
   spread_final(:, ivar, ilevel) = err_spread;

end
end % End of level loop

clear truth ens spread err err_spread

%----------------------------------------------------------------------
% Each variable in its own figure window
%----------------------------------------------------------------------
for ivar=1:pinfo.num_state_vars,

   figure(ivar); clf;
      ft       = netcdf(pinfo.truth_file);
      varunits = ft{pinfo.vars{ivar}}.units(:);
      close(ft);

      s1 = sprintf('%s %s Ensemble Mean for %s', model,pinfo.vars{ivar},pinfo.diagn_file);

      switch lower(pinfo.vars{ivar})
         case {'ps'}
            h1 = plot(times,          rms(:, ivar, 1), '-'); hold on;
            h2 = plot(times, spread_final(:, ivar, 1), '--');

            s{1} = sprintf('time-mean Ensemble Mean error  = %f', mean(         rms(:, ivar, 1)));
            s{2} = sprintf('time-mean Ensemble Spread = %f',      mean(spread_final(:, ivar, 1)));
         otherwise
            h1 = plot(times, squeeze(         rms(:, ivar, :)),'-'); hold on;
            h2 = plot(times, squeeze(spread_final(:, ivar, :)),'--');

            for i = 1:num_levels,
               s{i           } = sprintf('level %d error  %.3f', i,mean(         rms(:, ivar, i)));
               s{i+num_levels} = sprintf('level %d spread %.3f', i,mean(spread_final(:, ivar, i)));
            end
      end

      %h = legend([h1 h2],s); legend(h,'boxoff')
      h = legend(s); legend(h,'boxoff')
      grid on;
      xlabel(sprintf('time (%s) %d timesteps',timeunits,num_times))
      ylabel(sprintf('global-area-weighted distance (%s)',varunits))
      title(s1,'interpreter','none','fontweight','bold')

end



%----------------------------------------------------------------------
% helper function
%----------------------------------------------------------------------
function slice = GetPS(fname,copyindex, tstart, tend);

% this gets a bit funky if tstart == tend; in spite of trying to
% pass 0 in for the squeeze argument, it seems to squeeze anyway,
% leaving the dims in the wrong slots.  copyindex will always
% squeeze out, so try to figure out if you get a 3d or 2d return
% back from getnc() and then deal with time explicitly.
% same goes below in GetLevel().
corner     = [tstart, copyindex, -1, -1];
endpnt     = [tend,   copyindex, -1, -1];
ted        = getnc(fname,'ps',corner,endpnt);
n = ndims(ted);
if (n == 2)
   [ny,nx] = size(ted);
   nt = 1;
else
   [nt,ny,nx] = size(ted);
end
slice      = reshape(ted,[nt ny*nx]);


function slice = GetLevel(fname,varstring,copyindex,ilevel,tstart,tend);
%
corner     = [tstart, copyindex, ilevel, -1, -1];
endpnt     = [tend,   copyindex, ilevel, -1, -1];
ted        = getnc(fname,varstring,corner,endpnt);
n = ndims(ted);
if (n == 2)
   [ny,nx] = size(ted);
   nt = 1;
else
   [nt,ny,nx] = size(ted);
end
slice      = reshape(ted,[nt ny*nx]);


%----------------------------------------------------------------------
% helper function
%----------------------------------------------------------------------
function wts = SphereWeights(lons, lats)
% SphereWeights creates weights based on area ...
%
% lons,lats must be 1D arrays (in degrees)

nlons = length(lons);
nlats = length(lats);

if ( prod(size(lons)) ~= nlons )
   disp('longitude array is of higher dimension than anticipated.')
   error('Must be a vector.')
end
if ( prod(size(lats)) ~= nlats )
   disp('latitude array is of higher dimension than anticipated.')
   error('Must be a vector.')
end

rads    = zeros(nlats,1);               % Ensure lats is a column vector,
rads(:) = pi*lats/180.0;                % and convert to radians.
wts     = cos( rads ) * ones(1,nlons);  % Results in a [nlat-x-nlon] matrix.
wts     = wts ./ sum(wts(:));           % Normalize to unity.
