%% DART:random_proof checks the integrity of the random number generation
%                    under usage patterns common to DART.
%
%  Run test_reseed.f90; postprocess_test_reseed.csh 

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

seedsep  = [2, 60, 3600, 21600, 43200, 86400]; 
nseeds   = length(seedsep);
nreseeds = 1000;
xax      = 1:nreseeds;

u1.mean = zeros(nreseeds,nseeds);
u1.std  = zeros(nreseeds,nseeds);
g1.mean = zeros(nreseeds,nseeds);
g1.std  = zeros(nreseeds,nseeds);
gM.mean = zeros(nreseeds,nseeds);
gM.std  = zeros(nreseeds,nseeds);

for i = 1:nseeds,

   u = load(sprintf('seeds_%05d/unif_1000'    ,seedsep(i)));
   g = load(sprintf('seeds_%05d/gauss_1000'   ,seedsep(i)));
   G = load(sprintf('seeds_%05d/gauss_1000000',seedsep(i)));

   u1.mean(:,i) = u(:,1);
   u1.std( :,i) = u(:,2);
   g1.mean(:,i) = g(:,1);
   g1.std( :,i) = g(:,2);
   gM.mean(:,i) = G(:,1);
   gM.std( :,i) = G(:,2);

end

%% uniform(0,1) length 1000

figure(1); clf

nbins = 50;
hm = zeros(nbins,nseeds); xm = zeros(nbins,nseeds);
hs = zeros(nbins,nseeds); xs = zeros(nbins,nseeds);
legstr = cell(1,nseeds);

for i = 1:nseeds,
   [hm(:,i), xm(:,i)] = hist(u1.mean(:,i),nbins);
   [hs(:,i), xs(:,i)] = hist(u1.std( :,i),nbins);
   
   legstr{i} = sprintf('date separation %d',seedsep(i));
end

subplot(2,1,1)
   plot(xm,hm,'-*')
   h = legend(legstr);
   legend(h,'boxoff')
   title({'distribution of uniform(0,1) for a length = 1000 series',...
             'state vector times separated by 2 seconds'})
   ylabel('1000 seeds')
   xlabel('mean')

subplot(2,1,2)
   plot(xs,hs,'-*')
   ylabel('1000 seeds')
   xlabel('standard deviation')
   
orient tall
print -dpdf uniform_1000.pdf

%% gaussian(0,1) length 1000

figure(2); clf

for i = 1:nseeds,
   [hm(:,i), xm(:,i)] = hist(g1.mean(:,i),nbins);
   [hs(:,i), xs(:,i)] = hist(g1.std( :,i),nbins);
end

subplot(2,1,1)
   plot(xm,hm,'-*')
   h = legend(legstr);
   legend(h,'boxoff')
   title({'distribution of gaussian(0,1) for a length = 1000 series', ...
             'state vector times separated by 2 seconds'})
   ylabel('1000 seeds')
   xlabel('mean')
      
subplot(2,1,2)
   plot(xs,hs,'-*')
   ylabel('1000 seeds')
   xlabel('standard deviation')    

orient tall
print -dpdf gaussian_1000.pdf

   
%% gaussian(0,1) length 1000000
figure(3); clf

for i = 1:nseeds,
   [hm(:,i), xm(:,i)] = hist(g1.mean(:,i),nbins);
   [hs(:,i), xs(:,i)] = hist(g1.std( :,i),nbins);
end

subplot(2,1,1)
   plot(xm,hm,'-*')
   h = legend(legstr);
   legend(h,'boxoff')
   title({'distribution of gaussian(0,1) for a length = 1,000,000 series', ...
             'state vector times separated by 2 seconds'})
   ylabel('1000 seeds')
   xlabel('mean')
   
subplot(2,1,2)
   plot(xs,hs,'-*')
   ylabel('1000 seeds')
   xlabel('standard deviation')
   
orient tall
print -dpdf gaussian_1000000.pdf

%% plot of the first 100 values - gaussian(0,1) distrubution

datmat = zeros(100,1000);
n = 2;

for i = 1:nseeds,
    
   figure(i); clf
   
   fname = sprintf('seeds_%05d/ran_test_r8_%05d_numbers.out',seedsep(i),seedsep(i));
   datmat(:) = load(fname);
   first = datmat(1,:);
   dupes = find(first == first(1));
   fprintf('There are %d first values identical to the first value.\n',length(dupes))
   
   sname = sprintf('seeds_%05d/seeds',seedsep(i));
   bob   = load(sname);
   seeds = bob(:,2); clear bob;
   
   disp('These are the seeds that generate the identical first values:')
   for j = 1:length(dupes),
       fprintf('  %d\n',seeds(dupes(j)))
   end
   
   subplot('position',[0.1 0.1 0.76 0.8])
      plot(datmat(1:n,:)','.')
      h = title({sprintf('First %d values for each reseeding - gaussian(0,1) distribution',n),fname});
      set(h,'Interpreter','none')
      xlabel('reseeding number')
      ylabel('gaussian(0,1)')
      axlims = axis;
      
   subplot('position',[0.87 0.1 0.1 0.8])
      [N,Y] = hist(datmat(:),nbins);
      plot(N,Y)
      axis([-Inf Inf axlims(3) axlims(4)])
      set(gca,'YAxisLocation','right')
      
end

