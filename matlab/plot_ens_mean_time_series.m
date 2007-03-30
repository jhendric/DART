% DART : Plots time series of ensemble mean and truth
%                                                                               
% plot_ens_mean_time_series    interactively queries for the needed information.
%              Since different models potentially need different pieces of 
%              information ... the model types are determined and additional 
%              user input may be queried.
%
% Ultimately, plot_ens_mean_time_series will be replaced by a GUI.
%
%
% All the heavy lifting is done by PlotEnsMeanTimeSeries.

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

if (exist('truth_file') ~= 1)
   truth_file = input('Input name of True State file; <cr> for True_State.nc\n','s');
   if isempty(truth_file)
      truth_file = 'True_State.nc';
   end
end

if (exist('diagn_file') ~=1)
   disp('Input name of prior or posterior diagnostics file;')
   diagn_file = input('<cr> for Prior_Diag.nc\n','s');
   if isempty(diagn_file)
      diagn_file = 'Prior_Diag.nc';
   end
end

vars  = CheckModel(diagn_file);   % also gets default values for this model.

if (exist(truth_file)==2)
   CheckModelCompatibility(truth_file, diagn_file)
   truth_file = truth_file;
else
   truth_file = [];
end

switch lower(vars.model)

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
	 'lorenz_04','forced_lorenz_96','ikeda'} 

      varid = SetVariableID(vars);      % queries for variable IDs if needed.
      pinfo = struct('truth_file', truth_file, ...
                     'diagn_file', diagn_file, ...
                     'var'       , varid.var, ...
                     'var_inds'  , varid.var_inds);

      disp(sprintf('Comparing %s and \n          %s', pinfo.truth_file, pinfo.diagn_file))
      disp(['Using State Variable IDs ', num2str(pinfo.var_inds)])

   case 'fms_bgrid'

      pinfo = GetBgridInfo(diagn_file, 'PlotEnsMeanTimeSeries');
      pinfo.truth_file = truth_file;   % since it has been verified to be compatible.
      pinfo.diagn_file = diagn_file;   % since it has been verified to be compatible.

   case 'cam'

      vars.truth_file     = truth_file;
      vars.diagn_file     = diagn_file;
      vars.prior_file     = [];
      vars.posterior_file = [];
      pinfo               = GetCamInfo(vars, 'PlotEnsMeanTimeSeries');
   %  pinfo.copyindices   = SetCopyID2(vars.diagn_file);
   %  pinfo.copies        = length(pinfo.copyindices);
      pinfo.truth_file    = truth_file;
      pinfo.diagn_file    = diagn_file;

   otherwise

      error(sprintf('model %s not implemented yet', vars.model))

end

pinfo % echo for posterity.

PlotEnsMeanTimeSeries( pinfo )
clear vars varid
