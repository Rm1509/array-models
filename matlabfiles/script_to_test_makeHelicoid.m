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
homefol='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\230217\Set1_RangeOfAR';PN = 3;paramset = [1, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]; naparam = 'AR'; % this was used for the aspect ratio (single pitch length)
homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\230221test';
homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\difftnumbersTest';PN = 2;paramset = [1,2,3,4,5,9,10,12, 15, 16]; naparam = 'no. particles';% numbers of particles in grid
PN = 2;paramset = [1,4,9,16,25,36,49,64,81,100]; naparam = 'no. particles';% new numbers of particles
homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\angular disorder';PN = 8;paramset = [0,1,5,10,20,45,70]; naparam = 'sig. ang. disorder';% angular disorder 230223
PN = 8;paramset = [15,20,25,30,35,40]; naparam = 'sig. ang. disorder';% new angular disorder sigma
homefol = 'C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\sizedisorder';PN = 9;paramset = [0,0.1,0.2,0.3,0.4,0.5]; naparam ='sig. size dis.';% sigma of size disorder
homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\locationdisorder';PN = 10;paramset = [0,10,20,30,50,80]; naparam = 'sig. loc. dis.';% sigma of location disorder
homefol='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\arraySpacing';PN = 11;paramset = [-50, -20, 0, 10,50,100]; naparam= 'spacer'; % spacer added to grid
homefol = 'C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\hexagonalArrayNumbers';PN = 2;paramset = [1,7,13,19,25,37,43,61]; naparam = 'no. particles'; % new numbers of particles
homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\manypitches';PN = 4;paramset = [1,2,3,4,5]; naparam = 'no. pitches'; % numbers of pitches MORE IS TOO MANY
homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\multilayers';PN = 11;paramset = [1,2,3,4]; naparam = 'no. multilayers'; % number of multilayers
homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\parabhelicAngDis';PN = 8; paramset= [0,1,5,10,20,45,70]; naparam = 'sig. ang. disorder';% sigma angular disorder
homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\parabhelicNumParticles'; paramset = [1,7,13,19]; naparam = 'part no.';
% homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\2layersAngDis'; paramset =[0,5,10,15,20,25,30,40]; naparam ='sig. ang. dis.';
% homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\2layersLocDis'; paramset =[0,5,10,20,30,50]; naparam ='sig. loc. dis.';
% homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\2layersUDDis'; paramset =[0,5,10,20,30,50]; naparam ='sig. xloc. dis.';
% homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\2layersSizeDis'; paramset =[0,0.05,0.1,0.2,0.3,0.4]; naparam ='sig. size dis.';
homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\CerintheCombo';
% homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\sizes';paramset = [60,80,100,120,140]; naparam = 'width';
% homefol = 'C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\multirods';%PN = 3;paramset = [1,2,3,4,5,9,10,12, 15, 16]; naparam = 'no. particles';% numbers of particles in grid
homefol = 'C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\singles'; paramset = [1,2,3,4,5,6,7];
% homefol='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\parabhelicNumParticles'; %paramset
homefol ='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\CerintheCombo\smallbox';paramset=ones(10);
homefol = 'C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\singleslowangle'; paramset = [1,2,3,4,5,6,7];

newhomef = makelumhomefol(homefol);
homefol=newhomef;
%%

numsh = 1;      % number of particles in square grid
rotsig =10;%70;     %10; % sigma of random rotation around y & z axes independently (deg)
resize1 = 0.3;  %0.2; % sigma of random particle resize uniform in all 3 dimensions
locsig = 0;%150;%150;   %20; % sigma of particle distribution deviation from grid in y-z 
xloc = 0;%150;%250;     % sigma of particle distribution deviation from grid in x 
smoono=5;
ar = 4;         %3;
PITCHNO = 1;
thk = [2,2,3];
rad = 100;      %[100,40];%
spacer1=0;      % this adds space between the elements of the array in the multiparticle
layers=1;
grid = 'hex';   % 'hex' or 'square'
shape='helicoidp';%'ring';% NB FOR RING GIVE TWO RADII % 'rod';%'helicoid';%
saveon = [0,1]; %[0,1]; %saves the multiparticle set only

%straightver
% numsh = 7;rotsig =0; resize1 = 0;locsig = 0; xloc = 0;smoono=5;ar = 4;PITCHNO = 1;thk = [2,2,3];rad = 100;spacer1=0; layers=1;grid = 'hex'; shape='helicoidp';saveon = [0,1];

% paramset = ones(1,6);
for j=1:5
for i= 1:size(paramset,2)

param=paramset(i);

[stuff, volrec] =        multiparticle(homefol, numsh, ar, PITCHNO, thk, saveon, rad, rotsig, resize1, locsig, xloc, spacer1, layers, shape, grid);
% [fv_combined, volrec] = multiparticle(homefol,numsh, ar, PITCHNO, thk,  saveon, rad, rotsig, resize1,locsig, xloc, spacer1, layers, shape, grid)

end
end