function bins = rank_hist(ens, verif)
%RANK_HIST: Computes a rank histogram given time series of ensemble and verification 

% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

% Determine how long the time series is
num_times = size(ens, 1);
ens_size = size(ens, 2);

% Need ens_size + 1 bins
bins(1:ens_size + 1) = 0.0;

% Loop through time series to get count for each bin
for i = 1:num_times
   count = 0;
   for j = 1:ens_size
      if verif(i) > ens(i, j)
         count = count + 1;
      end
   end
   bins(count+1) = bins(count+1) + 1;
end



