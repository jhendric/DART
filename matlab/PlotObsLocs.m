% PLOTOBSLOCS
%
% Plot an observation_locations.NNN.dat file as output from the 
% latest obs_diag program in DART.  (You must enable a namelist entry
% to get the output files; they are not created by default.)
%
% warning: this has a very long argument list, because it is intended that 
% it be called from another program or function.   all the arguments can be 
% replaced by the character string 'default' to use the internal defaults.
%
% usage:
% PlotObsLocs(in_used, in_box, in_typelist, in_epochlist, in_subset, 
%             in_plotd, in_world, in_invertz, in_writeplot, 
%             in_legend2dloc, in_legend3dloc, in_viewlist)
%
% where:
%
% in_used =  0 to plot both used and unused observations
%         = -1 for unused obs only
%         =  1 for used obs only  (default)
% 
% in_box = optional bounding box [xmin xmax ymin ymax] for 2d or 
%          [xmin xmax ymin ymax zmin zmax] for 3d.  default is max extent of 
%          the data, or 0:360/-90:90 if in_world is true
% 
% in_typelist = list of integer observation types to plot.  all types
%               are plotted by default.  (hint: examine the end of 
%               the ObsDiagAtts.m file which is created by the obs_diag
%               program to see what number is used for each observation type.)
% 
% in_epochlist = list of integer time epochs to plot.  output from obs_diag 
%                program creates multiple files starting at 1 for each time 
%                period.  default is all.
% 
% in_subset = number of points to plot for each observation type.  a subset 
%             of the observations is randomly selected to thin down very 
%             large observation files.  the default is to plot all.
% 
% in_plotd = 2 or 3 for 2d or 3d plot, respectively.  default is 2d.
% 
% in_world = 0 for no, 1 for yes.  if true, plot the outlines of the continents
%            of the world on a lon/lat grid.  default is yes.
% 
% in_invertz = 0 for no, 1 for yes. if true, invert the Z axis.  pressure 
%              levels, for example, are largest at the earth's surface and 
%              decrease as you go up.  default is no.
% 
% in_writeplot = 0 for no, 1 for yes.   if true, write out each plot in a 
%                color postscript file.  default is no.  warning -- this can
%                be slow and the files can be large for very large numbers of
%                observations.
% 
% in_legend2dloc = where to put the legend box which shows the plot markers 
%                  for each observation type.  the default values are often 
%                  directly over a region of interest and so far it has seemed
%                  impossible to find a completely satisfactory default.  
%                  type 'help legend' to get a list of possible strings.  
%                  this is used for a 2d plot.
% 
% in_legend3dloc = where to put the legend box which shows the plot markers 
%                  for each observation type.  the default values are often 
%                  directly over a region of interest and so far it has seemed
%                  impossible to find a completely satisfactory default.  
%                  type 'help legend' to get a list of possible strings.  
%                  this is used for a 3d plot.
% 
% in_viewlist = [az el; az el; az el; ...] list of viewpoints for 3d plotting.
%               the default is 2 views of each plot at [10  10;  10  80 ]. 
%               [0 90] [90 0] and [0 0] may be interesting options.
%               the 3d plots can be rotated interactively and the [az el] is 
%               printed in the lower left corner during the rotation; these
%               views can be used to create .ps files if you are saving the 
%               plots to a file or if you want the same viewpoints reused.
% 
% 

%------------------------------------------------------
function PlotObsLocs(in_used, in_box, in_typelist, in_epochlist, in_subset, in_plotd, in_world, in_invertz, in_writeplot, in_legend2dloc, in_legend3dloc, in_viewlist)

% data subset selections:
%  used vs unused (or both) obs
%  subregions - lon, lat, height min/max
%  obs type list
%  random subset to thin data
%  epoch list
%
% plot type and output format:
%  2d or 3d plot  
%  not lat/lon (no world map)
%  create .ps files for output or not


% defaults for everything that has code to handle alternatives below.
% suggestion - leave the defaults alone and in the next section overwrite
% the values with alternative values.

arg_used = 1;         % -1 = unused, 0 = both, 1 = used
arg_box = [];         % [[lon_min, lon_max, lat_min, lat_max], v_min, v_max]
arg_typelist = [];    % numeric observation type list; if ~[], integer list
arg_epochlist = [];   % list of epochs to process; if ~[], integer list
arg_subset = 0;       % random subset to thin data; if > 0, how many to keep
arg_plotd = 2;        % 2 = 2d, 3 = 3d
arg_world = 1;        % worldwide (lon/lat) locations, plot over map (0=no)
arg_invertz = 0;      % for 3d plot, 1=z axis plots N to 0, 0 = 0 to N
arg_writeplot = 0;    % write out postscript plot files (1=yes)
arg_legend2dloc = 'SouthWest';  % where on the plot to plop the 2d legend box
arg_legend3dloc = 'NorthWest';  % where on the plot to plop the 3d legend box
                              % see matlab doc for alternatives
arg_viewlist = [10  10;  10  80; ];
                      % default viewpoints for 3D plots.
                      % [90 0], [0 90] and [0 0] are also good.



% get the values from the arguments and fill them in:
if (~isa(in_used,'char'))
  arg_used = in_used;
end
if (~isa(in_box,'char'))
  arg_box = in_box;
end
if (~isa(in_typelist,'char'))
  arg_typelist = in_typelist;
end
if (~isa(in_epochlist,'char'))
  arg_epochlist = in_epochlist;
end
if (~isa(in_subset,'char'))
  arg_subset = in_subset;
end
if (~isa(in_plotd,'char'))
  arg_plotd = in_plotd;
end
if (~isa(in_world,'char'))
  arg_world = in_world;
end
if (~isa(in_invertz,'char'))
  arg_invertz = in_invertz;
end
if (~isa(in_writeplot,'char'))
  arg_writeplot = in_writeplot;
end
if (in_legend2dloc ~= 'default')
  arg_legend2dloc = in_legend2dloc;
end
if (in_legend3dloc ~= 'default')
  arg_legend3dloc = in_legend3dloc;
end
if (~isa(in_viewlist,'char'))
  arg_viewlist = in_viewlist;
end


%-----------


% if user specified an epoch list, use only those.
% else, set it to a huge number and bail when we
% run out of files to open.
if (~isequal(arg_epochlist,[]))
  epochlist = arg_epochlist;
else
  epochlist = 1 : 100000;
end


% for each epoch file given or that we find:
for epoch = epochlist

 % get data in and ignore the 1 header line.  using textread brings in
 % the data as a single long 1-D array, so i use reshape and then transpose 
 % to get it into the original [N 6] matrix shape.
 
 % this assumes the data is in a series of files named 
 % 'observation_locations.NNN.dat', where NNN is epoch 001, 002, etc
 % the latest updates to the obs_diag program create these files
 % if the 'print_obs_locations' namelist entry is set to .true.
 %
 
 datafile = sprintf('observation_locations.%03d.dat', epoch);
 
 try
  r = textread(datafile, '%f', 'headerlines', 1);
 catch
  % last file found
  break
 end
 
 % format of the data in these files must be:
 %   lon  lat  ivert  flavor  key  used
 t = reshape(r, 6, []);
 s = transpose(t);
 
 % start each plot clean
 clf; 

 % load the strings generated at the same time as the locations
 % this includes the text description of what the observation was
 % in the 'Observation_Kind()' array.  (this file is produced by
 % running obs_diag -- full name is ObsDiagAtts.m)
 
 ObsDiagAtts;
 
 % up to 36 different observation types, using a different marker
 % and color for each.
 markers = { 
  'r.', 'b.', 'g.', 'c.', 'm.', 'y.', 
  'r+', 'b+', 'g+', 'c+', 'm+', 'y+', 
  'r*', 'b*', 'g*', 'c*', 'm*', 'y*', 
  'ro', 'bo', 'go', 'co', 'mo', 'yo', 
  'rx', 'bx', 'gx', 'cx', 'mx', 'yx', 
  'rd', 'bd', 'gd', 'cd', 'md', 'yd', 
 };
 
 
 % setup before looping over observations:
 mobs = max(s(:,4));       % max obs type number found in file
 nobstypes = 0;            % running count of obs types found
 lmax = 1;                 % largest level found (default 1 so not degenerate)
 obs_labels = {};          % legend labels
 
 % if we are only looking at specific observation types loop over those.
 % otherwise generate a list from 1 to the largest kind number found in 
 % the input file and loop over those.
 if (~isequal(arg_typelist,[]))
   obslist = arg_typelist;
 else
   obslist = 1 : mobs;
 end
 
 % loop over observation types, plotting each in 2d or 3d with
 % a different marker.
 for obs = obslist
 
   % select only this obs type, and cycle if there are none
   thisobs = s(s(:, 4) == obs, :);
   if (isequal(thisobs,[]))
     continue;
   end
 
   % select the rows we want to plot, optionally with a variety
   % of subsetting.  function form is:
   %   out = select_subset(rawdata_in, used/unused, boundingbox, subset_count)
   l = select_subset(thisobs, arg_used, arg_box, arg_subset);
 
   % if there are any of this obs type left, plot them
   if (size(l, 1) ~= 0) 
     if (arg_plotd == 2) 
       % 2d plot of (lon,lat):
       phandle = plot(l(:,1),l(:,2),markers{obs}); 
     else
       % 3D plot of (lon,lat,level):
       phandle = plot3(l(:,1),l(:,2),l(:,3),markers{obs});
       thismax = max(l(:,3));
       lmax = max(lmax, thismax);    % save overall max height for axis
     end
 
     hold on;   % accumulate on same plot
 
     % add this observation type to the legend string array
     nobstypes = nobstypes + 1;
     obs_labels(nobstypes) = Observation_Kind(obs) ;
 
   end    % if data of this observation type exists
 
 end  %  for obslist
 
 
 % if this is a world map plot, add in an outline map of
 % the continents and set default size to 0/360,-90/90 (which can 
 % be overridden with the box argument).  if not world, default to
 % the data min/max or the box argument.
 
 % use the first 4 items for a 2d plot, all (6) for 3d;
 % arg_box is the user override; use_box is either the
 % users choice or our default.
 use_box = [];
 if (~isequal(arg_box, []))
    if (arg_plotd == 2)
       use_box = arg_box(1:4);
    else
       if (size(arg_box,2) > 4)
          use_box = arg_box;
       else
          use_box = [arg_box 0 lmax];
       end
       if (arg_invertz == 1)
           set(gca,'ZDir','reverse');
       end
    end
 else
    if (arg_world)
      if (arg_plotd == 2)
        use_box = [0 360 -90 90];
      else
        use_box = [0 360 -90 90 0 lmax];
        if (arg_invertz == 1)
           set(gca,'ZDir','reverse');
           %set(gca,'ZScale','log');
        end
      end
    end
 end
 
 % get a handle to the axes graphics object
 % ahandle = axes

 % if we actually set something, use it to constrain the axis limits.
 if (~isequal(use_box, []))
     axis(use_box);
 else
     axis tight;
 end
    
 % set legend, and try to shrink the original size of the legend bounding box
 % because it is pretty large by default.   a 'good' location depends on the
 % kind of plot, and even so this may be something you want to override
 % depending on where in the plot the interesting data falls.
 if (arg_plotd == 2)
   legendloc = arg_legend2dloc;
 else
   legendloc = arg_legend3dloc;
 end

 h = legend( obs_labels , 'Location', legendloc, 'Interpreter', 'none', 'FontSize', 8);
 
 % example of how to escape only underscores if we still want to use tex
 % in the strings. (instead of turning the interpreter off completely).
 %h=title(strrep(Observation_Kind(obs),'_','\_'));
  
 % whole world in (lon, lat) degree coords
 if (arg_world)
 
    % add a 2D plot of the world continent outlines
lmax
    worldmap(lmax);
    
    % these plots are generally longer than high, and add 3d-box.
    orient landscape;
    set(gca, 'Box', 'on');

    % various attempts to make the x/y axis have the same spacing per degree
    % but shrink the vertical because it can be much larger.  but all these
    % did strange things to the viewpoint so i am giving up on them for now.
    %axis equal;

    %vert = lmax / 100;
    %set(gca, 'DataAspectRatio', [1 1 vert], 'PlotBoxAspectRatio', [1 1 vert], 'Box', 'on');
    %set(gca, 'DataAspectRatio', [1 1 10], 'Box', 'on');
 
   xlabel('Longitude (degrees)', 'FontSize', 14);
   ylabel('Latitude (degrees)', 'FontSize', 14);
   if (arg_plotd == 3)
     zlabel('Height (units = pressure?)', 'FontSize', 14);
   end
 else
   xlabel('First coordinate', 'FontSize', 14);
   ylabel('Second coordinate', 'FontSize', 14);
   if (arg_plotd == 3)
     zlabel('Third coordinate', 'FontSize', 14);
   end
 end
    
 % add input filename somewhere?
 % text()

 % make it look roughly like it will when printed.
 wysiwyg;
 
 
 if (arg_used < 0)
   tstring = sprintf('Unused Observation Locations at Epoch %d', epoch);
 elseif (arg_used > 0)
   tstring = sprintf('Used Observation Locations at Epoch %d', epoch);
 else
   tstring = sprintf('Used + Unused Observation Locations at Epoch %d', epoch);
 end

 if (arg_subset > 0)
    nstring = sprintf(', %d random obs per type', arg_subset);
    tstring = strcat(tstring, nstring);
 end

 title(tstring, 'FontSize', 16);
 
 
 % if 3d plot, turn on grid and view for interactive look
 if (arg_plotd == 3)
   grid on;
   % set up viewlist as: arg_viewlist = [ az el; az el; ... ];
   for v = 1 : size(arg_viewlist, 1)
      view(arg_viewlist(v, :));

      disp('Pausing.  Hit any key to continue');
      pause;
   end
 else
   disp('Pausing.  Hit any key to continue');
   pause;
 end
 
 % save plots into files.  2d has only 1 plot; 3d has many possible
 % view angles -- you can set a list and get multiple plots out.
 if (arg_writeplot) 
    if (arg_plotd == 2)
      fname = sprintf('2d_locations_epoch%.03d.ps', epoch);
      print('-dpsc', fname);
    else

      % set up viewlist as: arg_viewlist = [ az el; az el; ... ];
      for v = 1 : size(arg_viewlist, 1)
         view(arg_viewlist(v, :));
         fname = sprintf('3d_locations%d_epoch%.03d.ps', v, epoch);
         print('-dpsc', fname);
      end

    end   % 2d vs 3d
 end  % create .ps files?
 
end  % for epoch  (main loop)

disp('Plots done.');

end  % function PlotObsLoc3D


%------------------------------------------------------
% given a raw data array and possible selectors, return
% the subset of rows which match.

function out = select_subset(raw, used, box, number)

% default is to assume we found no data.  set it at the
% end if we get there.
out = [];

% make sure some data is passed in before trying to subset it
if (isequal(raw, []))
  fprintf('no observations passed into select_subset\n');
  return
end


% select used, unused, or all observations (column 6)
if (used < 0)
  data = raw(raw(:,6)<1, :);
elseif (used > 0)
  data = raw(raw(:,6)>0, :);
else
  data = raw;
end

if (isequal(data, []))
  fprintf('all observations removed after used/unused selection\n');
  return
end

% select subset by area if one is given
if (~isequal(box, []))
  data = data(data(:, 1) > box(1), :);  % lon min/max
  data = data(data(:, 1) < box(2), :);
  data = data(data(:, 2) > box(3), :);  % lat min/max
  data = data(data(:, 2) < box(4), :);
  if (size(box,2) > 4)
    data = data(data(:, 3) > box(5), :);  % vertical if given
    data = data(data(:, 3) < box(6), :);
  end
end

if (isequal(data, []))
  fprintf('all observations removed after bounding box selection\n');
  return
end


% select a random subset of the remaining values 
numleft = size(data, 1);
if ((number > 0) && (numleft > number))
  % but cannot have repeats in list
  sel = sort(pickme(number, numleft));
  if (~isequal(sel, []))
    data = data(sel, :);
  end
end

if (isequal(data, []))
  fprintf('all observations removed after random subset selection\n');
  return
end

out = data;

end % function out = select_subset(raw, used, box, number)

%------------------------------------------------------
% generate a list of R random numbers between 1 and N
% (inclusive) with *no* repeats.  Return them in sorted order.
function list = pickme(R, N)
 
% assume failure until assured otherwise.
list = [];

% test for infinite loop 
if (R > N) 
  fprintf('cannot select %d unique values between 1 and %d\n', R, N);
  return
end

% ok, get started.
selections = zeros(R, 1);
picked = zeros(N, 1);
next = 0;  
 
while (next < R)
  % generate a random candidate, but if it has already been
  % selected, cycle.
  candidate = ceil(N .* rand);
  if (picked(candidate) == 1)
    continue;
  end

  % running count of how many we have
  next = next + 1;
  selections(next) = candidate;
  picked(candidate) = 1;

end

list = selections;

end % function list = pickme(R, N)

