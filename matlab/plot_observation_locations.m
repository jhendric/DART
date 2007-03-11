% PLOT_OBSERVATION_LOCATIONS : Plots the locations of the input observations
%
% By default this command creates 2d plots of observation locations,
% one per time epoch, from data output from the obs_diag program if
% the 'print_obs_locations' namelist item in the &obs_diag list is .true.
%
% There are lots of user settable options.  This script prompts you
% interactively for the most common ones.  Then it calls PlotObsLocs()
% with the proper argument list to pass in your selections.
%
% See the documentation for PlotObsLocs() -- it has a lot of arguments in the
% calling sequence.

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

% setup all args to be the string 'default', which will be interpreted by 
% the PlotObsLocs routine to use the default values.   
 
plotd       = 'default';
used        = 'default';
typelist    = 'default';
box         = 'default';
epochs      = 'default';
subset      = 'default';
world       = 'default';
writeplot   = 'default';
loc2dstring = 'default';
loc3dstring = 'default';
viewlist    = 'default';
invertz     = 'default';
 
% what the arg list looks like:
%PlotObsLocs(in_used, in_box, in_typelist, in_epochlist, in_subset, in_plotd, in_world, in_invertz, in_writeplot, in_legend2dloc, in_legend_3dloc, in_viewlist)

done = 0;
disp('Plot observations at their proper locations.  Many subsetting options exist.');
disp('Hitting <cr> to answer the questions will use the default value,');
disp('or once you have made a selection, reuse the previous value.');
disp(' '); 
disp('The default plotting options are:');
disp('  2D plot, full world map, all obs types, all times, ');
disp('  no file output, Z axis increases up.');
disp(' '); 

% loop and keep the previous default until the user says to quit
while done == 0

   % 2D or 3D plot?
   reply = input('Input 2 for 2D plot, 3 for 3D plot:  ');
   if (~isempty(reply))
      plotd = reply;
   end
    
   % plot used, unused, or both
   reply = input('Input -1=unused obs, 0=both, 1=used:  ');
   if (~isempty(reply))
      used = reply;
   end
   
   % restrict observations to a particular observation type?
   reply = input('Input [obs type list] to plot only some obs types, ''default'' to reset:  ');
   if (~isempty(reply))
      typelist = reply;
   end
    
   % restrict observations to a particular subregion?
   reply = input('Input [xmin xmax ymin ymax] for bounding box, ''default'' to reset:  ');
   if (~isempty(reply))
      box = reply;
   end
    
   % restrict input to particular time epochs?
   reply = input('Input [epoch list] for particular times, ''default'' to reset:  ');
   if (~isempty(reply))
      epochs = reply;
   end
    
   % sample data to reduce counts?
   reply = input('Input count for random subset of each obs type, ''default'' to reset:  ');
   if (~isempty(reply))
      subset = reply;
   end
    
   % plot world map beneath?
   reply = input('Input 0 to remove world map, 1 to restore it:  ');
   if (~isempty(reply))
      world = reply;
   end
    
   % write out .ps files?
   reply = input('Input 1 to write .ps files for each plot, 0 to reset:  ');
   if (~isempty(reply))
      writeplot = reply;
   end
    
   % legendloc
   reply = input('Input Matlab string for legend location, ''default'' to reset:  ');
   if (~isempty(reply))
      if (plotd == 3)
         loc3dstring = reply;
      else
         loc2dstring = reply;
      end
   end
    
   % viewlist and invert z axis
   invertz = 1;
   if (plotd == 3) 
      reply = input('Input [az el; az el] list for 3d views, ''default'' to reset:  ');
      if (~isempty(reply))
         viewlist = reply;
      end
   
      reply = input('Input 1 to invert Z axis (e.g. for pressure); 0 otherwise:  ');
      if (~isempty(reply))
         invertz = reply;
      end
   end
    
   PlotObsLocs(used, box, typelist, epochs, subset, plotd, world, invertz, writeplot, loc2dstring, loc3dstring, viewlist)
   
   reply = input('Quit now or plot again? 1 = quit, <cr> plots again:  ');
   if (~isempty(reply))
      done = reply;
   end
  
end

