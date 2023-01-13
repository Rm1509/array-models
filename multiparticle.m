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


% this is required if you want to save it,  but is not actually important
homefol = 'C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\230112';

%this calls the basic shape, for when the same shape is reproduced multiple
%times
eghandle = makeHelicoid(3,  1,   [2,2,3], homefol, 0,  100, 0.05);

%figure setup
figure
pcol = parula(9);

% to make a multiarray, each one is made in this loop
for i = 1: 9
    
    %place a shape in space at a random location
    p = patch('Vertices', eghandle.vertices+repmat([150,200*randn(1,2)], size(eghandle.vertices,1),1), 'Faces', eghandle.faces);
    hold on

    % indicate the colour of the shape
    p.FaceColor= pcol(i,:);
    p.EdgeColor= pcol(i,:);

    %rotate the shape according to a distribution 
    rotate(p, [randn(1), randn(1), randn(1)], 45.*randn(1))
    
end

% neaten up the plot
looks(gcf)


end

function looks(ffig)
    ffig    
    view(3); 
    axis tight
    camlight 
    lighting gouraud
end