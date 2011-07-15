function x_new = advance_oned(x, alpha, model_bias)
%% advance_oned(x, alpha, model_bias)

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

x_new = x + comp_dt(x, alpha, model_bias);
end

%---------------------------------------------------

% Internal function comp_dt
function dx = comp_dt(x, alpha, model_bias)

% Compute the time tendency; alpha controls nonlinearity
% model_bias controls a shift in the model dynamics
dx = (x + model_bias) + alpha .* x .* abs(x);

end
