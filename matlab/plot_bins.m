%% DART:plot_bins Plots ensemble rank histograms 
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

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% error_checking ... 
% exist('bob') == 1   means the variable exists. 
%                     the value of the variable is checked later.

if (exist('truth_file','var') ~= 1) 
   truth_file = input('Input name of True State file; <cr> for True_State.nc\n','s');
   if isempty(truth_file)
      truth_file = 'True_State.nc';
   end
end

if (exist('diagn_file','var') ~=1)
   disp('Input name of prior or posterior diagnostics file;') 
   diagn_file = input('<cr> for Prior_Diag.nc\n','s');
   if isempty(diagn_file)
      diagn_file = 'Prior_Diag.nc';
   end
end

if ( exist(truth_file,'file') ~= 2 ), error('%s does not exist.',truth_file); end
if ( exist(diagn_file,'file') ~= 2 ), error('%s does not exist.',diagn_file); end

pinfo = CheckModelCompatibility(truth_file,diagn_file);
vars  = CheckModel(truth_file);   % also gets default values for this model.

switch lower(vars.model)

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
	 'lorenz_04', 'forced_lorenz_96','ikeda','simple_advection'}

      varid = SetVariableID(vars);      % queries for variable IDs
      pinfo.truth_file = truth_file;
      pinfo.diagn_file = diagn_file;
      pinfo.var        = varid.var;
      pinfo.var_inds   = varid.var_inds;
      clear varid

   case 'fms_bgrid'

      pinfo = GetBgridInfo(pinfo, diagn_file, 'PlotBins');

   case 'cam'

      pinfo = CombineStructs(pinfo,vars);
      pinfo = GetCamInfo(pinfo, diagn_file, 'PlotBins');

   case 'pe2lyr'

      pinfo = GetPe2lyrInfo(pinfo, diagn_file, 'PlotBins');

   case 'mitgcm_ocean'

      pinfo = GetMITgcm_oceanInfo(pinfo, diagn_file, 'PlotBins');

   otherwise

      error('model %s not implemented yet', vars.model)

end

PlotBins(pinfo);
clear vars
