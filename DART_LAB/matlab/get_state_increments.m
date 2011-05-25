function [state_incs] = get_state_increments(state_ens, obs_ens, obs_incs)
%% get_state_increments Computes state increments given observation increments and
% the state and obs prior ensembles

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% Compute state variance and covariance
covar = cov(state_ens, obs_ens);


state_incs = obs_incs * covar(1, 2) / covar(2, 2);

