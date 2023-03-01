% script to test makeHelicoid

% this is one place that homefol has been sent 'C:\Users\Rox\OneDrive -
% University of Bristol\Documents\lumerical beb\221201 multi tests\matfiles'

%previously
% eghandle = makeHelicoid(3,  1,      [2,2,3],  1, homefol, 100, 0.05  );
% eghandle = makeHelicoid(3,  1,   [2,2,3], homefol, 1,  100, 0.05);
% eghandle = makeHelicoid(3,  1,   [2,2,3], homefol, 0,  100, 0.05);


%now

multiparticle(newhomef)
%then run lumerical file with same
%%
homefol='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\230217\Set1_RangeOfAR';
homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\230221test';
homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\difftnumbersTest';
homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\angular disorder';
homefol = 'C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\sizedisorder';
homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\locationdisorder';
homefol='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\arraySpacing';
homefol = 'C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\hexagonalArrayNumbers';
homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\manypitches';

newhomef = makelumhomefol(homefol);
homefol=newhomef;
%%
paramset = [1, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]; % this was used for the aspect ratio (single pitch length)
paramset = [1,2,3,4,5,9,10,12, 15, 16]; % this was used for numbers of particles in grid
paramset = [0,0.1,0.2,0.3,0.4,0.5,0.75,1];
paramset = [0,1,5,10,20,45,70]; % this was the set for angular disorder 230223
paramset = [1,4,9,16,25,36,49,64,81,100]; % new numbers of particles
paramset = [15,20,25,30,35,40]; % new angular disorder digma
paramset = [0,0.1,0.2,0.3,0.4,0.5]; % sigma of size disorder
paramset = [0,10,20,30,50,80]; % sigma of location disorder
paramset = [-50, -20, 0, 10,50,100]; % spacer added to grid
paramset = [1,7,13,19,25,37,43,61]; % new numbers of particles
paramset = [1,2,4,6,8,10];

% paramset =1;
numsh = 19;  % number of particles in square grid
rotsig =0;  %10; % sigma of random rotation around y & z axes independently (deg)
resize1 = 0;%0.2; % sigma of random particle resize uniform in all 3 dimensions
locsig = 0; %20; % sigma of particle distribution deviation from grid
smoono=5;
ar = 4;%3;
PITCHNO = 10;
thk = [2,2,3];
rad = 100;
spacer1=0;  % this adds space between the elements of the array in the multiparticle

saveon = [0,1];%[0,1]; %saves the multiparticle set only

for i= 1:size(paramset,2)

param=paramset(i);

[stuff, volrec] =        multiparticle(homefol, numsh, ar, param, thk, saveon, rad, rotsig, resize1, locsig,spacer1);
% [fv_combined, volrec] = multiparticle(homefol,numsh, ar, PITCHNO, thk,  saveon, rad, rotsig, resize1,locsig, spacer1)
end
