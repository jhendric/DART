function PlotCorrel(fname, base_var_index, base_time)
% Plots space-time series of correlation between a given variable at a given
% time and all other variable at all times in an ensemble time sequence.

if (exist(fname) ~= 2), error(sprintf('%s does not exist.',fname)), end

% Get some file-specific information.
f = netcdf(fname,'nowrite');
model      = f.model(:);
num_vars   = ncsize(f{'StateVariable'}); % determine # of state variables
num_times  = ncsize(f{'time'}); % determine # of output times
num_copies = ncsize(f{'copy'}); % determine # of ensemble members
close(f)

% disp(sprintf('PlotCorrel: fname is %s',fname))
% disp(sprintf('PlotCorrel: base_var_index is %d',base_var_index))
% disp(sprintf('PlotCorrel: base_time      is %d',base_time))
% disp(sprintf('PlotCorrel: num_vars is %d',num_vars))

% The Base Variable Index must be a valid state variable

if ( base_var_index > num_vars )
   disp( sprintf('%s only has %d state variables', fname, num_vars))
   error(sprintf('you wanted variable # %d ', base_var_index))
end

% The Time must be within range also.

if ( base_time > num_times )
   disp( sprintf('%s only has %d output times', fname, num_times))
   error(sprintf('you wanted time # %d ', base_time))
end

statevariables = getnc(fname,'StateVariable');

% Get 'standard' ensemble series 
base = get_ens_series(fname, base_var_index);

% It is efficient to preallocate correl storage ... 
correl = zeros(num_vars,num_times);

% Need to loop through all variables in the ensemble
for i = 1:num_vars
   state_var = get_ens_series(fname, i);
   correl(i, :) = ens_correl(base, base_time, state_var);
end

% Now for the plotting part ...
clf;

contour(correl,[-1:0.2:1]);
s1 = sprintf('%s Correlation of state variable %d, T = %d of %s', ...
         model, base_var_index, base_time,fname);
s2 = sprintf('against all variables, all times, all %d ensemble members', ...
         num_copies-2); 
title({s1,s2},'interpreter','none','fontweight','bold')
xlabel('time (timestep #)')
ylabel('state variable (index)')
set(gca,'YTick',statevariables)
colorbar

% highlight the reference state variable and time

hold on;
plot(base_time,base_var_index,'kh','MarkerSize',12,'MarkerFaceColor','k')

% Y = { '\tt 1', '\tt 2', '\bf 3', '\tt 4', '\tt 5', '\tt 6', '\tt 7', '\tt 8', '\tt 9'}; set(gca,'YTickLabel',Y);
