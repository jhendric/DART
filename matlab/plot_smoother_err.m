% DART : summary plots of global error and spread using the smoother
% Example 1
% diagn_file = 'Posterior_Diag.nc';
% truth_file = 'True_State.nc';   % for smoother, is Lag_00001_Diag.nc better?
% num_lags   = 10;
% plot_total_err

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

lag_file   = 'Lag_%05d_Diag.nc'; % pattern for lag file names

if (exist('num_lags')   ~= 1), num_lags = 10000; end
if (exist('truth_file') == 1), 
   def_true = truth_file;
else
   def_true = 'True_State.nc';
end

disp('Input name of True State file;')
truth_file = input(sprintf('<cr> for %s\n',def_true),'s');
if isempty(truth_file)
     truth_file = def_true;
end

% Loop over all possible lags, if the corresponding netCDF file 
% does not exist, we automatically terminate.

for lag=1:num_lags

  def_diag = sprintf(lag_file, lag);
  
  disp('Input name of smoother lag diagnostics file;')
  diagn_file = input(sprintf('<cr> for %s\n', def_diag),'s');
  if isempty(diagn_file)
     diagn_file = def_diag;
  end

  if (exist(diagn_file) ~= 2)
     disp('file does not exist. Must be done.')
     return
  end

  pinfo = CheckModel(diagn_file);
  pinfo.truth_file = truth_file;
  pinfo.diagn_file = diagn_file;

  bob = CheckModelCompatibility(truth_file, diagn_file);
  pinfo.truth_time = bob.truth_time;
  pinfo.diagn_time = bob.diagn_time;

  clear bob

  disp(sprintf('Comparing %s and \n          %s', ...
                pinfo.truth_file, pinfo.diagn_file))
  
  PlotTotalErr( pinfo );

  disp(' ')

end

clear pinfo
