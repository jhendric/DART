function hop_test_check(file0, file1, file2)
% DART hop_test_check 
%
% USAGE:
% x = hop_test_check(file0, file1, file2);
%
% file0    is the filename of the initial state 
% file1    is the filename of the 'long hop' case
% file2    is the filename of the 'short hop' case
%
% The difference between file2 and file1 will be scaled by the 
% time tendency - the amount the field changed between file0 and file1.
%
% Example:
%
% file0 = '/gpfs/ptmp/thoar/restarts/CAM/caminput_1.nc';
% file1 = '/glade/scratch/thoar/archive/hop_24/rest/2008-11-01-00000/hop_24.cam_0001.i.2008-11-01-00000.nc';
% file2 = '/glade/scratch/thoar/archive/hop_12/rest/2008-11-01-00000/hop_12.cam_0001.i.2008-11-01-00000.nc';
% x     = hop_test_check(file0, file1, file2);

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if ( (exist(file0,'file') ~= 2) || ...
     (exist(file1,'file') ~= 2) || ...
     (exist(file2,'file') ~= 2) )
   fprintf('One of the input files does not exist.\n')
   fprintf('%s, \n',file0)
   fprintf('%s, or\n',file1)
   error('%s',file2)
end

%% extract all the variable names from both files and find the intersection.
vars1 = sort(FindProgVars(file1));
vars2 = sort(FindProgVars(file2));
vars  = sort(FindCommonVars(vars1,vars2));

%% loop over the variables and extract them from both files.
%  take the difference and convert it to a percentage of 'file1-file0' ... 
% Since Matlab automatically squeezes out singleton dimensions, and 'time'
% is usually a singleton in these files, I need to squeeze the singleton
% dimension out of the variable dimension strings, too.

bob = jet128; colormap(bob);

for i = 1:length(vars)

   fprintf('Comparing %s\n',vars{i})
    
   varinfo       = nc_getvarinfo(file0,vars{i});
   nonsingletons = (varinfo.Size > 1);
   mydimnames    = varinfo.Dimension(nonsingletons);
   mydimsizes    = varinfo.Size(nonsingletons);
   levdim        = find(strcmp(mydimnames,'lev'));
   
   start         = nc_varget(file0,vars{i});
   onehop        = nc_varget(file1,vars{i});
   twohop        = nc_varget(file2,vars{i});
   tendency      = onehop - start;
   change        = twohop - onehop;

   % Need to also plot vars that do not have a 'lev' dimension
   % PS only defined on surface, for example.
   
   for ilevel = 1:mydimsizes(levdim)
       myplot(vars{i}, ilevel, change, tendency)
       fprintf('          level %d ... \n',ilevel)
       pause(0.1)
   end
   
end


function myplot(varname,levelindx,diffmat,tendmat)
%% Make some plots
% To be comparable with the other figures in the article, the hemispheres 
% need to be [-180,180], the major ticks and labels are every 30 degrees.
% no titles.
% bigger fonts
% Peter wanted .eps or .ps - good - can use 'painters'

slab = squeeze(diffmat(levelindx,:,:));
slabmin = min(slab(:));
slabmax = max(slab(:));
clim = [slabmin slabmax];

orgslab = squeeze(tendmat(levelindx,:,:));
orgmin = min(orgslab(:));
orgmax = max(orgslab(:));
datmax = max(abs([orgmin orgmax]));
clim = [-datmax datmax];

if slabmin == slabmax, return; end

figure(1); clf; orient tall; 
subplot(2,1,1)
   imagesc(slab,clim);
   set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
   title(sprintf('%s level %d difference from hopping',varname,levelindx))
   axis image
   h = colorbar('vert');
   
   subplot(2,1,2)
   imagesc(orgslab,clim);
   set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
   title(sprintf('%s level %d tendency',varname,levelindx))
   axis image
   h = colorbar('vert');






function vars = FindProgVars(fname)
%% Determine the multi-dimensional variables that are not coordinate variables.
% These are probably the variables of interest.

fileinfo  = nc_info(fname);
nvars     = length(fileinfo.Dataset);
isPROGvar = zeros(nvars,1);

% We are not interested in the 1D variables, so just skip them.
for i = 1:nvars
   if ( length(fileinfo.Dataset(i).Size) > 1 ), isPROGvar(i) = 1; end
end

if (sum(isPROGvar) == 0)
   error('No multidimensional variables in %s',fname)
end

% coerce just the names into a cell array 

varind = 0;
for i = 1:nvars
   if (isPROGvar(i) > 0)
      varind = varind + 1;
      vars{varind} = fileinfo.Dataset(i).Name;
   end
end



function vars = FindCommonVars(vars1,vars2)
% 
k = 0;

for i = 1:length(vars1)
   if ( any(strcmp(vars1{i},vars2)) )
      k = k + 1;
      vars{k} = vars1{i};
   end
end


function x = jet128()
x = jet(128);
x(64:65,:) = 1;


