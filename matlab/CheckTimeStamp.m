  function CheckTimeStamp(directory, varargin)
% print the timestamps of BINARY DART-format files. 
%
% Example 1: (check all instances of filter_ics.nnnn)
%
% CheckTimeStamp('/glade/scratch/syha/work_hires')
% 
% Example 2: (binary file of the wrong 'endian'-ness)
% directory = '/glade/scratch/syha/work_hires';
% endian = 'little';
% CheckTimeStamp(directory,'endian','little')


%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% Helen, aim: parse inputs, bulletproofing
% Inputs: 
%   Required
%      directory - directory contating .ics files 
%   Optional:
%      endian - ieee-be (big) or ieee-le (liitle) 
%      fbase - base of ics files, e.g. filter_ics
% outputs: prints to screen filename and timestamp

p = inputParser; % create new instance of parser class
p.FunctionName = 'input parser :: requires input directory (string); valid optional inputs are, endian (string), fbase (string)';

%% set defaults for optional parameters
defaultEndian   = 'native';
defaultFbase    = 'filter_ics';

addRequired(  p, 'directory', @ischar); % require directory, check input is a character array
addParamValue(p,    'endian', defaultEndian,   @ischar);
addParamValue(p,     'fbase', defaultFbase,    @ischar);

p.parse(directory, varargin{:}) % parse inputs

%% collect the results of parsing (makes code easier to read)

endianIn = p.Results.endian;
fbase    = p.Results.fbase;

%% check inputs

assert(exist(directory, 'dir')==7, 'directory %s does not exist', directory)
   
tempstr  = strcat(fbase, '*');
ens_size = length(dir(fullfile(directory, tempstr)));

if (ens_size > 0)
   fprintf('The ensemble size is %d\n',ens_size)
else
   error('no %s files exist in %s',fbase,directory);
end

switch lower(endianIn)
   case {'big','ieee-be'} 
      endian = 'ieee-be';
   case {'little','ieee-le'} 
      endian = 'ieee-le';
   otherwise
      endian = 'native';
end
   
timebase = datenum(1601,1,1); % start of Gregorian Calender

for i = 1:ens_size

   fname = sprintf('%s/%s.%04d',directory,fbase,i);
   
   if exist(fname, 'file')
   
       fid   = fopen(fname,'rb',endian); 
       rec1  = fread(fid,4,'int32');

       days = rec1(3);
       seconds = rec1(2);
       fdays = seconds/86400; 
       datestring = datestr(timebase + days + fdays);

       fprintf('%s.%04d has timestamp of %d %d    %s\n',fbase,i,days,seconds,datestring)
       fclose(fid);
   
   else
       fprintf('WARNING : %s.%04d does not exist \n',fbase,i)
   end

end

