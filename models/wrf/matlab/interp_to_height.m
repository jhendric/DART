function var_interp = interp_to_height( var_in, heights, level)
% interp_to_height Interpolates to a height level given heights
%
% var_interp = interp_to_height( var_in, heights, level)
%
% Interpolates var_in(:,:,:) to a height level given heights
% in heights(:,:,:).  Interpolation is linear in height. 
% Set var_interp to NaN where level is beneath (<) heights(1,:,:).

% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

 [Nk Nj Ni] = size(heights) ; 

 below = zeros(Nj,Ni) ;  var_below = zeros(Nj,Ni) ;
 dz =    zeros(Nj,Ni) ;       dvar = zeros(Nj,Ni) ;

 % at each horiz. location, find highest height beneath level
 for ii = 1:Ni; for jj = 1:Nj; 
    % kk(jj,ii) = max( find( heights(:,jj,ii) - level < 0 ) ) ;
    kk = max( find( heights(:,jj,ii) - level < 0 ) ) ;
    % if isempty(kk); disp([ kk , jj, ii ] ); end
    if isempty(kk); kk = 1; end  % level is below surface here

    below(jj,ii) = heights(max(1 ,kk  ),jj,ii) ;  % height below level
    dz(jj,ii) = heights(min(Nk,kk+1),jj,ii)    ...
                 - heights(max(1 ,kk  ),jj,ii) ;
    var_below(jj,ii) = var_in(max(1 ,kk  ),jj,ii) ;
    dvar (jj,ii) = var_in(min(Nk,kk+1),jj,ii) - var_in(max(1 ,kk  ),jj,ii) ;
 end; end

 var_interp =  ( dvar ./ dz ) .* ( level - below )  + var_below ;

 var_interp( level < heights(1,:,:) ) = NaN ;
   % level is beneath lowest mass level
