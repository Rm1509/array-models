function multiparticle()



% fff = figure; 
homefol = 'C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\230112';
eghandle = makeHelicoid(3,  1,   [2,2,3], homefol, 0,  100, 0.05);

% p=patch('Vertices', eghandle.vertices, 'Faces', eghandle.faces);

% p.FaceColor= 'cyan';
% p.EdgeColor= 'cyan';
figure
pcol = parula(9);
for i = 1: 9

p = patch('Vertices', eghandle.vertices+repmat([150,200*randn(1,2)], size(eghandle.vertices,1),1), 'Faces', eghandle.faces);
p.FaceColor= pcol(i,:);
p.EdgeColor= pcol(i,:);
% rotate(p, [1,1,1].*randn(1,3), 45.*randn(1))
rotate(p, [randn(1), randn(1), randn(1)], 45.*randn(1))
hold on
end
% aspect([1 1 1])
looks(gcf)


end
function looks(ffig)
ffig    
view(3); 
    axis tight
    camlight 
    lighting gouraud
end