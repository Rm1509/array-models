% script to test makeHelicoid

% this is one place that homefol has been sent 'C:\Users\Rox\OneDrive -
% University of Bristol\Documents\lumerical beb\221201 multi tests\matfiles'

%previously
% eghandle = makeHelicoid(3,  1,      [2,2,3],  1, homefol, 100, 0.05  );
% eghandle = makeHelicoid(3,  1,   [2,2,3], homefol, 1,  100, 0.05);
% eghandle = makeHelicoid(3,  1,   [2,2,3], homefol, 0,  100, 0.05);


%now

% multiparticle(newhomef)
%then run lumerical file with same
%%
PN = 2;paramset = [1,4,9,16,25,36,49,64,81,100]; naparam = 'no. particles';% new numbers of particles
PN = 8;paramset = [15,20,25,30,35,40]; naparam = 'sig. ang. disorder';% new angular disorder sigma

% This file uses relative paths, so will save outputs in a 'shapes_data' directory
% in the project root
abspath = fileparts(mfilename('fullpath')); % Can't run as a live document
homefol = fullfile(abspath, '../shapes_data');
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
