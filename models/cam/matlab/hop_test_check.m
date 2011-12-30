function hop_test_check(file0, file1, file2, varname)
% DART hop_test_check 
%
% USAGE:
% hop_test_check(file0, file1, file2);
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
% hop_test_check(file0, file1, file2)
%
% Example:
%
% file3 = '/gpfs/ptmp/thoar/restarts/CLM/clminput_1.nc';
% file4 = '/glade/scratch/thoar/archive/hop_24/rest/2008-11-01-00000/hop_24.clm2_0001.r.2008-11-01-00000.nc';
% file5 = '/glade/scratch/thoar/archive/hop_12/rest/2008-11-01-00000/hop_12.clm2_0001.r.2008-11-01-00000.nc';
% hop_test_check(file3, file4, file5)

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

if (nargin > 3)
    vars = {varname};
    pausecmd = 'fprintf(''pausing at level %d ... hit a key to continue\n'',ilevel), pause';
else
    pausecmd = 'fprintf(''           level %d ... \n'',ilevel); pause(0.1)';
end

for i = 1:length(vars)
    
   varinfo = nc_getvarinfo(file0,vars{i});
   if (varinfo.Nctype == 2)
       % Character string variables need not be checked.
       fprintf('Skipping   %s\n',vars{i})
       continue
   else
       fprintf('Comparing  %s\n',vars{i})
   end
   nonsingletons = (varinfo.Size > 1);
   mydimnames    = varinfo.Dimension(nonsingletons);
   mydimsizes    = varinfo.Size(nonsingletons);
   levdim        = find(strcmp(mydimnames,'lev'));

   myunits       = GetAttribute(file0,vars{i},'units');
   start         = nc_varget(file0,vars{i});
   onehop        = nc_varget(file1,vars{i});
   twohop        = nc_varget(file2,vars{i});
   tendency      = onehop - start;
   change        = twohop - onehop;
   
   if (length(mydimsizes) == 1)
       continue % How do we compare 1D arrays?
   elseif (isempty(levdim))
       % We still have a multimensional object... 
       my2dplot(vars{i}, change, tendency, mydimnames, myunits)
       disp('          pausing - hit any key to continue ...')
       pause
   else
       for ilevel = 1:mydimsizes(levdim)
           myplot(vars{i}, ilevel, change, tendency, mydimnames, myunits)
           eval(pausecmd)
       end  
   end
 
end


function my2dplot(varname,slab,orgslab,dimnames,units)
%% Make some plots
%

slabmin = min(slab(:));
slabmax = max(slab(:));

orgmin = min(orgslab(:));
orgmax = max(orgslab(:));
datmax = max(abs([orgmin orgmax]));

if orgmin == orgmax
    clim = [-1 1];
else
    clim = [-datmax datmax];
end

sbpos = [0.10 0.06 0.80 0.28; 
         0.10 0.43 0.80 0.28;
         0.10 0.80 0.80 0.15]; 

% Need to know what dimensions we are plotting

figure(1); clf; orient tall; 

subplot('position',sbpos(3,:))
   hist(slab(:),50)
   title(sprintf('(min %0.5g %s) hopping difference histogram (max %0.5g %s)',slabmin,units, slabmax,units))
   xlabel(units)

subplot('position',sbpos(2,:))
   imagesc(slab,clim);
   set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
   title(sprintf('%s difference from hopping',varname))
   axis image
   ylabel(sprintf('%s (index)',dimnames{1}))
   h = colorbar('vert');
   set(get(h,'YLabel'),'String',units)
   
subplot('position',sbpos(1,:))
   imagesc(orgslab,clim);
   set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
   title(sprintf('%s tendency',varname))
   axis image
   ylabel(sprintf('%s (index)',dimnames{1}))
   xlabel(sprintf('%s (index)',dimnames{2}))
   h = colorbar('vert');
   set(get(h,'YLabel'),'String',units)


function myplot(varname,levelindx,diffmat,tendmat,dimnames,units)
%% Make some plots
%

slab = squeeze(diffmat(levelindx,:,:));
slabmin = min(slab(:));
slabmax = max(slab(:));

if slabmin == slabmax, return; end

orgslab = squeeze(tendmat(levelindx,:,:));
orgmin = min(orgslab(:));
orgmax = max(orgslab(:));
datmax = max(abs([orgmin orgmax]));
clim = [-datmax datmax];

sbpos = [0.10 0.06 0.80 0.28; 
         0.10 0.43 0.80 0.28;
         0.10 0.80 0.80 0.15]; 

figure(1); clf; orient tall; 
subplot('position',sbpos(3,:))
   hist(slab(:),50)
   title(sprintf('(min %0.5g %s) hopping difference histogram (max %0.5g %s)',slabmin,units,slabmax,units))
   xlabel(units)

subplot('position',sbpos(2,:))
   imagesc(slab,clim);
   set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
   title(sprintf('%s level %d difference from hopping',varname,levelindx))
   axis image
   ylabel(sprintf('%s (index)',dimnames{2}))
   h = colorbar('vert');
   set(get(h,'YLabel'),'String',units)
   
subplot('position',sbpos(1,:))
   imagesc(orgslab,clim);
   set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
   title(sprintf('%s level %d tendency',varname,levelindx))
   axis image
   ylabel(sprintf('%s (index)',dimnames{2}))
   xlabel(sprintf('%s (index)',dimnames{3}))
   h = colorbar('vert');
   set(get(h,'YLabel'),'String',units)



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

function attvalue = GetAttribute(fname,varname,attname)
varinfo = nc_getvarinfo(fname,varname);
attvalue = [];
for i = 1:length(varinfo.Attribute)
   if (strcmp(varinfo.Attribute(i).Name,attname))
      attvalue = varinfo.Attribute(i).Value;
      break
   end
end
