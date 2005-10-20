% DART : Plots ensemble rank histograms 
%
% plot_bins    interactively queries for the information needed to create
%              ensemble rank histograms. Since different models potentially
%              need different pieces of information ... the model types are
%              determined and additional user input may be queried.
%
% Ultimately, plot_bins will be replaced by a GUI.
% In the end, the heavy lifting is done by PlotBins.
%
% Example 1 (for low-order models)
%
% truth_file = 'True_State.nc';
% diagn_file = 'Prior_Diag.nc';
% plot_bins

% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$
 
% error_checking ... 
% exist('bob') == 1   means the variable exists. 
%                     the value of the variable is checked later.

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

CheckModelCompatibility(truth_file,diagn_file)
vars  = CheckModel(truth_file);   % also gets default values for this model.
varid = SetVariableID(vars);      % queries for variable IDs if needed.

switch lower(vars.model)

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
	 'lorenz_04', 'forced_lorenz_96'}

      pinfo = struct( 'truth_file'    , truth_file, ...
                      'diagn_file'    , diagn_file, ...
                      'var'           , varid.var , ...
                      'state_var_inds', varid.var_inds);

   %  disp(sprintf('Comparing %s and \n          %s', pinfo.truth_file, pinfo.diagn_file))
   %  disp(sprintf('Using Variable %s IDs %s', pinfo.var,num2str(pinfo.state_var_inds)))

   case 'fms_bgrid'

      pinfo = GetBgridInfo(diagn_file, 'PlotBins');
      pinfo.truth_file = truth_file;   % since it has been verified to be compatible.
      pinfo.diagn_file = diagn_file;   % since it has been verified to be compatible.

   otherwise

      error(sprintf('model %s not implemented yet', vars.model))

end

pinfo

PlotBins(pinfo);
clear vars varid
