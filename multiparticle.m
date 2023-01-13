function multiparticle()
% 
% -	Neaten up multiparticle so all of these things are inputs and outputs
% -	Allow or disallow intersection
% -	Control distribution in z – first attempt at this, but it doesn’t work very well
% -	Add a grid of distribution in x-y which the deviation is away from
% -	Control the angular orientation within a range
% -	Save these and load them in Lumerical (NB. Radically different sizes!)
% -	Control size and add these to Lumerical
% -	Allow a range of sizes and a range of lengths
% -	Twist the curls round the central corkscrew
% -	Use different shapes

%how many shapes wanted
numsh = 9;

%figure setup
figure
pcol = parula(numsh);

% this is required if you want to save it,  but is not actually important
homefol = 'C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\230112';

%this calls the basic shape, for when the same shape is reproduced multiple
%times
eghandle = makeHelicoid(3,  1,   [2,2,3], homefol, 0,  100, 0.05);

%how big is the object
obsz = obsz3(eghandle)

%make an array
arlen = ceil(sqrt(numsh));
arym = repmat([1:arlen]', 1, arlen)*obsz(2,3);
arzm = repmat([1:arlen], arlen, 1)*obsz(3,3);
arxm = zeros(size(arym))*obsz(1,3)/10;
arxyz = [arxm(:), arym(:), arzm(:)];
figure; scatter3(arxm(:), arym(:), arzm(:))

% to make a multiarray, each one is made in this loop
for i = 1: numsh
    
    %choose spatial location for point
    arxyzi = arxyz(i,:);
    edloc = [-obsz(1,1), 0,0];
    arloc = arxyzi;
    varloc = [0,randn(1,2)]*0;
%     size(arxyzi)
%     size(edloc)
%     size(arloc)
%     size(varloc)

    newloc = repmat(arloc + edloc+ varloc, size(eghandle.vertices,1),1);
%     newloc = repmat([150,200*randn(1,2)], size(eghandle.vertices,1),1); % this is random

    
    %place a shape in space at a random location
    p = patch('Vertices', eghandle.vertices+newloc, 'Faces', eghandle.faces);
    hold on

    % indicate the colour of the shape
    p.FaceColor= pcol(i,:);
    p.EdgeColor= pcol(i,:);

    %rotate the shape according to a distribution 
%     rotate(p, [randn(1), randn(1), randn(1)], 45.*randn(1))
    
end

% neaten up the plot
looks(gcf)


end

%function to find size of 3d object
function obsz = obsz3(eghandle)
    for i = 1:3
        [jx,j] = min(eghandle.vertices(:,i));
        [kx,k] = max(eghandle.vertices(:,i));
        obsz(i,1:2)= [jx,kx];
        obsz(i,3) = kx-jx;
        % hold on; scatter3(eghandle.vertices(j,1),eghandle.vertices(j,2),eghandle.vertices(j,3),'r')
    end
end

%function to neaten plot
function looks(ffig)
    ffig    
    view(3); 
    axis tight
    camlight 
    lighting gouraud
    xlabel('x')
    ylabel('y')
    zlabel('z')
end