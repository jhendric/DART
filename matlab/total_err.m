function err = total_err(pred, verif, w)
%TOTAL_ERR: Computes Total error for time series of set of state variables

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

% Pred and verif are time_series_length x number of variables
num_times = size(pred, 1);
num_vars  = size(pred, 2);

if ( nargin == 2 )
   %--------------------------------------------------
   % unweighted error
   %--------------------------------------------------

   err = sqrt( sum( (pred - verif).^2       , 2) );

elseif ( nargin == 3 )
   %--------------------------------------------------
   % weighted error
   %--------------------------------------------------

   if (length(w) ~= num_vars)
      error(sprintf('wrong sizes %d and %d',length(w),num_vars))
   end

   err = zeros(num_times,1);

   for i = 1:num_times,

      d  = (pred(i,:) - verif(i,:)).^2;

      err(i) = sqrt( sum( d(:) .* w ) );
   end

else
   error('Wrong number of arguments.')
end