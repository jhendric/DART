% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$
 
fname = 'Prior_Diag';
tlon = getnc(fname, 'XLON');
we = size(tlon, 2);
tlat = getnc(fname, 'XLAT');
sn = size(tlat, 1);
level = getnc(fname, 'level');
bt = size(level, 1);
copy = getnc(fname, 'copy');
ens_size = size(copy, 1) - 2;

true_times = getnc(fname, 'time');
num_true_times = size(true_times, 1);
%num_true_times = 16

mean_ind = ens_size + 1;
std_ind = mean_ind + 1;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)])

for field_num = [1:9]

start_var = 1;
nx = we + 1;
ny = sn;
var_units = 'U (m/s)';
maxlev = bt;
maxval = 8.0;
if field_num > 1
start_var = start_var + bt*(we + 1)*sn ;
nx = we;
ny = sn + 1;
var_units = 'V (m/s)';
end
if field_num > 2
start_var = start_var + bt*we*(sn + 1);
nx = we;
ny = sn;
var_units = 'W (m/s)';
maxlev = bt + 1;
maxval = 0.03;
end
if field_num > 3
start_var = start_var + (bt + 1)*we*sn;
var_units = 'GZ (m^2/s^2)';
maxval = 700.0;
end
if field_num > 4
start_var = start_var + (bt + 1)*we*sn;
var_units = 'T (K)';
maxlev = bt;
maxval = 4.0;
end
if field_num > 5
start_var = start_var + bt*we*sn;
var_units = 'MU (Pa)';
maxlev = 1;
maxval = 700.0;
end
if field_num > 6
start_var = start_var + we*sn;
var_units = 'QV (kg/kg)';
maxlev = bt;
maxval = 0.001;
end
if field_num > 7
start_var = start_var + bt*we*sn*(field_num-7);
var_units = 'QC (kg/kg)';
maxval = 0.00003;
end
if field_num > 8
var_units = 'QR (kg/kg)';
maxval = 0.00007;
end

x = [1:2*num_true_times];
rmse = x;
spread = x;
E2 = x-x;

rms_mem = zeros(2*num_true_times,ens_size);

f_size = nx*ny*maxlev;

end_var = start_var + f_size - 1;

for itime = 1:num_true_times

% Extract field

fname = 'True_State';
field_vec_truth = getnc(fname, 'state',[itime -1 start_var],[itime -1 end_var],[1 1 1]);

fname = 'Prior_Diag';
field_vec_prior = getnc(fname, 'state',[itime mean_ind start_var],[itime mean_ind end_var],[1 1 1]);

field_vec = getnc(fname, 'state',[itime std_ind start_var],[itime std_ind end_var],[1 1 1]);

spread(2*itime-1) = sqrt((field_vec'*field_vec)/(f_size));

field_vec = field_vec_prior - field_vec_truth;
rmse(2*itime-1) = sqrt((field_vec'*field_vec)/(f_size));

x(2*itime-1) = true_times(itime,1)*24; % time in hours

fname = 'Posterior_Diag';
field_vec_posterior = getnc(fname, 'state',[itime mean_ind start_var],[itime mean_ind end_var],[1 1 1]);

field_vec = getnc(fname, 'state',[itime std_ind start_var],[itime std_ind end_var],[1 1 1]);

spread(2*itime) = sqrt((field_vec'*field_vec)/(f_size));

field_vec = field_vec_posterior - field_vec_truth;
rmse(2*itime) = sqrt((field_vec'*field_vec)/(f_size));

x(2*itime) = true_times(itime,1)*24; % time in hours

  for imem = 1:ens_size

fname = 'Prior_Diag';
field_vec_prior = getnc(fname, 'state',[itime imem start_var],[itime imem end_var],[1 1 1]);
field_vec = field_vec_prior - field_vec_truth;
rms_mem(2*itime-1,imem) = sqrt((field_vec'*field_vec)/(f_size));
    
fname = 'Posterior_Diag';
field_vec_posterior = getnc(fname, 'state',[itime imem start_var],[itime imem end_var],[1 1 1]);
field_vec = field_vec_posterior - field_vec_truth;
rms_mem(2*itime,imem) = sqrt((field_vec'*field_vec)/(f_size));

E2(2*itime-1) = E2(2*itime-1) + rms_mem(2*itime-1,imem);
E2(2*itime) = E2(2*itime) + rms_mem(2*itime,imem);

  end

end


subplot(3,3,field_num);

%if (ens_size > 0.0)
%     E2 = E2/ens_size;
%     plot(x,rmse,x,spread,x,E2,'--m')
%else
     plot(x,rmse,x,spread)
%end

     plot_title = [var_units];

     title(plot_title)

%     xlabel('hours')

     room = (x(2*num_true_times)-x(1))/10;

%     axis ([(x(1)-room) (x(2*num_true_times)+room) min(min(rmse,spread)) max(max(rmse,spread))])
     axis ([(x(1)-room) (x(2*num_true_times)+room) 0.0 maxval])

end

%if (ens_size > 0.0)
%     legend('RMS error','Spread','E2')
%else
     legend('RMS error','Spread')
%end
xlabel('hours')

% Loop for another try
%map_wrf;

