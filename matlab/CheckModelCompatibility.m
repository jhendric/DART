function pinfo_out = CheckModelCompatibility(arg1, arg2)
%% CheckModelCompatibility tries to ensure that two netcdf files can be compared.
% There are 2 ways to call this:  with 2 filenames, or with an already existing
% pinfo struct (with 2 filenames and 2 2-vector arrays for start/stop times).
% This routine fills in the 2-vectors with the time overlap region in a
% pinfo struct.
% If the time indices are common between the 2 files it returns the 
% [start,count] indices for each array (indexing starts at 1,N).
% It is an error situation if there is no overlap ([-1,-1] for both). 

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$


if (nargin == 1)      % better be a pinfo struct with at least these fields
  file1 = arg1.truth_file;  % string
  file2 = arg1.diagn_file;  % string
  pinfo_out = arg1;
elseif (nargin == 2)  % pair of filenames
  file1 = arg1;             % truth_file
  file2 = arg2;             % diagn_file
  pinfo_out.truth_file = file1;
  pinfo_out.diagn_file = file2;
else
  error('Wrong number of arguments: must be 1 (pinfo) or 2 (file1,file2)')
end

if ( exist(file1,'file') ~= 2 ), error('(file1) %s does not exist.',file1); end
if ( exist(file2,'file') ~= 2 ), error('(file2) %s does not exist.',file2); end

% set this up for later
pinfo_out.truth_time = [-1,-1];
pinfo_out.diagn_time = [-1,-1];

%% Get some information from the file1
tmodel  = nc_attget(file1,nc_global,'model');

if (isempty(tmodel)) 
   error('%s has no ''model'' global attribute.',file1)
end

tnum_copies = dim_length(file1,'copy');
tnum_times  = dim_length(file1,'time');
times       = nc_varget( file1,'time');
timeunits   = nc_attget( file1,'time','units');
timebase    = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin  = datenum(timebase(1),timebase(2),timebase(3));
ttimes      = times + timeorigin;

[tnum_vars,tdims] = ModelDimension(file1,tmodel);
if ( tnum_vars <= 0 )
   error('Unable to determine resolution of %s.',file1)
end

%% Get some information from the file2
dmodel  = nc_attget(file1,nc_global,'model');

if (isempty(dmodel)) 
   error('%s has no ''model'' global attribute.',file2)
end

dnum_copies = dim_length(file2,'copy');
dnum_times  = dim_length(file2,'time');
times       = nc_varget( file2,'time');
timeunits   = nc_attget( file2,'time','units');
timebase    = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin  = datenum(timebase(1),timebase(2),timebase(3));
dtimes      = times + timeorigin;

[dnum_vars,ddims] = ModelDimension(file2,dmodel);
if ( dnum_vars <= 0 )
   error('Unable to determine resolution of %s.',file2)
end

% rudimentary bulletproofing
if (strcmp(tmodel,dmodel) ~= 1)
   fprintf('%s has model %s\n',file1,tmodel)
   fprintf('%s has model %s\n',file2,dmodel)
   error('no No NO ... models must be the same')
end
pinfo_out.model = tmodel;

if (prod(tnum_vars) ~= prod(dnum_vars))
   fprintf('%s has %d state variables\n',file1,prod(tnum_vars))
   fprintf('%s has %d state variables\n',file2,prod(dnum_vars))
   error('no No NO ... both files must have same shape of state variables.')
end

% if the lengths of the time arrays did not match, this used to be an
% error.  now we call a function to try to find any overlapping regions
% in the time arrays and pass them back up to the called in the pinfo struct.
% they then get used to extract the corresponding hyperslabs of data for
% the matching times.

% construct the pinfo struct in this function
pinfo_out = timearray_intersect(pinfo_out, file1, file2, ttimes, dtimes);

% fail here if the times had nothing in common.
if ( ( pinfo_out.truth_time(1) == -1 ) || ...
     ( pinfo_out.truth_time(2) == -1 ) || ...
     ( pinfo_out.diagn_time(1) == -1 ) || ...
     ( pinfo_out.diagn_time(2) == -1 ))
   fprintf('%s has %d timesteps, from %f to %f\n', ...
                file1,tnum_times,ttimes(1), ttimes(tnum_times))
   fprintf('%s has %d timesteps, from %f to %f\n', ...
                file2,dnum_times,dtimes(1), dtimes(dnum_times))
   error('These files have no timesteps in common')
end


%----------------
% min1,max1 and min2,max2 are the index numbers of the intersection of the
% two input arrays.  -1s in those numbers means no intersection.  1, length()
% means identical (could add a separate flag to simplify the calling code).
function pret = timearray_intersect(pinfo, file1, file2, times1, times2)

% for floating point comparisons, must be within this (single precision)
% roundoff
epsilon = 0.0000001;

% default return; no intersection
pret = pinfo;
pret.truth_file = file1;
pret.diagn_file = file2;
pret.truth_time = [-1,-1];
pret.diagn_time = [-1,-1];

% ensure times are increasing and monotonic, and do they need to be
% a constant delta or not?  compute delta array and validate those match?
% (to within an epsilon with floating pt roundoff)

% check for the no-brainer case - identical time arrays.
% watch out for the floating point compares, and the min/max are probably
% redundant with the (1) and (l) comparisons, but until we put in checks
% for monotonicity, it's a cheap safety check.
len = length(times1);
if (   (length(times1) == length(times2)) ...
    && (abs(min(times1) - min(times2)) < epsilon) ...
    && (abs(max(times1) - max(times2)) < epsilon) ...
    && (times1(1) == times2(1)) ...
    && (times1(len) == times2(len)))
  pret.truth_time = [1,len];   % start/count
  pret.diagn_time = [1,len];   % start/count
  return
end

% A is whichever array has the lower min.  this reduces the number of
% cases below we have to check for.
if (min(times1) < min(times2))
  A = times1;
  B = times2;
else
  B = times1;
  A = times2;
end

% precompute the data max, min, lengths using the A,B assignments
% also, if differences are < epsilon, force equality to simplify
% the comparison code below
lenA = length(A);
lenB = length(B);
minA = min(A);
minB = min(B);
maxA = max(A);
maxB = max(B);
if (abs(minA - minB) < epsilon) , minB = minA; end
if (abs(maxA - minB) < epsilon) , minB = maxA; end
if (abs(maxA - maxB) < epsilon) , maxB = maxA; end

% case 1: disjoint regions; simply return here because 
% return struct was initialized to the 'no intersection' case.
if ((minA < minB) && (maxA < minB))  
  return
end

% case 2: B fully contained in A; return corresponding index nums of overlap
% include equal start & end points in this case.
if ((minA <= minB) && (maxB <= maxA)) 
  minI = find(abs(A - minB) < epsilon);
  maxI = find(abs(A - maxB) < epsilon);
  minJ = 1;
  maxJ = lenB;
else
% case 3: partial overlap, minA lower than minB
  minI = find(abs(A - minB) < epsilon);
  maxI = lenA;
  minJ = 1;
  maxJ = find(abs(B - maxA) < epsilon);
end

% now map back to the original input order - this test must match exactly
% the one used initially to assign A and B above.
if (min(times1) < min(times2))
  min1 = minI;
  max1 = maxI;
  min2 = minJ;
  max2 = maxJ;
else
  min1 = minJ;
  max1 = maxJ;
  min2 = minI;
  max2 = maxI;
end

% now put the indices in the return struct and we are done.
pret.truth_time = [min1, max1-min1+1];   % start,count
pret.diagn_time = [min2, min2-max2+1];   % start,count

% return here


function x = dim_length(fname,dimname)

y = nc_isvar(fname,dimname);
if (y < 1) 
   error('%s has no %s dimension/coordinate variable',fname,dimname)
end
bob = nc_getdiminfo(fname,dimname);
x   = bob.Length;



function [x,y] = ModelDimension(fname,modelname)
% Check the base geometry of the grid
x = 0;
y = NaN;

switch lower(modelname)

   case 'wrf'
      diminfo = nc_getdiminfo(fname,  'west_east_d01'); dnum_lons = diminfo.Length;
      diminfo = nc_getdiminfo(fname,'south_north_d01'); dnum_lats = diminfo.Length;
      diminfo = nc_getdiminfo(fname, 'bottom_top_d01'); dnum_lvls = diminfo.Length;
         x = 3;
         y = [dnum_lons dnum_lats dnum_lvls];

   case 'cam'
      dnum_lons = dim_length(fname,'lon');
      dnum_lats = dim_length(fname,'lat');
      dnum_lvls = dim_length(fname,'lev');
         x = 3;
         y = [dnum_lons dnum_lats dnum_lvls];

   case 'pe2lyr'
      dnum_lons = dim_length(fname,'lon');
      dnum_lats = dim_length(fname,'lat');
      dnum_lvls = dim_length(fname,'lev');
         x = 3;
         y = [dnum_lons dnum_lats dnum_lvls];

   case 'fms_bgrid'
      dnum_lons = dim_length(fname,'TmpI');
      dnum_lats = dim_length(fname,'TmpJ');
      dnum_lvls = dim_length(fname,'lev');
         x = 3;
         y = [dnum_lons dnum_lats dnum_lvls];

   case 'mitgcm_ocean'
      dnum_lons = dim_length(fname,'XG');
      dnum_lats = dim_length(fname,'YG');
      dnum_lvls = dim_length(fname,'ZG');
         x = 3;
         y = [dnum_lons dnum_lats dnum_lvls];

   case 'lorenz_96_2scale'
      dnum_X = dim_length(fname,'Xdim');
      dnum_Y = dim_length(fname,'Ydim');
	 x = 2;
         y = [dnum_X dnum_Y];

   case 'simple_advection'
      y = dim_length(fname,'loc1d');
	 x = 1;

   otherwise
      y = dim_length(fname,'StateVariable');
	 x = 1;

end

