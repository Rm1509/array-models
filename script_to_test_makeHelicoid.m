% script to test makeHelicoid

% this is one place that homefol has been sent 'C:\Users\Rox\OneDrive -
% University of Bristol\Documents\lumerical beb\221201 multi tests\matfiles'

%previously
% eghandle = makeHelicoid(3,  1,      [2,2,3],  1, homefol, 100, 0.05  );
% eghandle = makeHelicoid(3,  1,   [2,2,3], homefol, 1,  100, 0.05);
% eghandle = makeHelicoid(3,  1,   [2,2,3], homefol, 0,  100, 0.05);


%now
homefol='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev';
makelumhomefol(homefol)
multiparticle(homefol)
%then run lumerical file with same