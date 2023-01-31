function fv_combined = multiparticle()
% 
% -	Control distribution in z – done
% -	Add a grid of distribution in x-y which the deviation is on - done 
% -	Neaten up multiparticle so all of these things are inputs and outputs
% -	Twist the curls round the central corkscrew
% -	Control the angular orientation within a range
% -	Allow a range of sizes and a range of lengths

% -	Allow or disallow intersection

% -	Save these and load them in Lumerical (NB. Radically different sizes!)
% -	Control size and add these to Lumerical

% -	Use different shapes

%how many shapes wanted
numsh = 2; % number of particles in square grid
rotsig =10; % sigma of random rotation around y & z axes independently (deg)
resize1 = 0.2; % sigma of random particle resize uniform in all 3 dimensions
locsig = 20; % sigma of particle distribution deviation from grid

ar = 3;
PITCHNO = 1;
thk = [2,2,3];
saveon = [0,1];
rad = 100;
einterval=0.05;
homefol = 'C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\230112';

%figure setup
% figure
pcol = parula(numsh);

% this is required if you want to save it,  but is not actually important

%this calls the basic shape, for when the same shape is reproduced multiple
%times
[f1,v1] = makeHelicoid(ar, PITCHNO, thk, homefol, saveon(1), rad, einterval);
% handle1 = makeHelicoid(3,  1,   [2,2,3], homefol, 0,  100, 0.05);
% ar=3; PITCHNO=1; thk = [2,2,3]; saveon=0;  rad = 100; einterval=0.05; homefol='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumericalbeb\221201 multi tests\matfiles'; 

close(gcf);
%how big is the object
obsz = obsz3(v1)

%make an array
arlen = ceil(sqrt(numsh));
arym = repmat([1:arlen]', 1, arlen)*obsz(2,3);
arzm = repmat([1:arlen], arlen, 1)*obsz(3,3);
arxm = zeros(size(arym))*obsz(1,3)/10;
arxyz = [arxm(:), arym(:), arzm(:)];
figure; scatter3(arxm(:), arym(:), arzm(:))

% fv_combined = struct([]);
fv_combined.vertices = [];
fv_combined.faces = [];
% to make a multiarray, each one is made in this loop

for i = 1: numsh

%     resizing
    xmul = 1 +randn(1)*resize1;%.*obsz(1,3);
    ymul = 1 +randn(1)*resize1;%.*obsz(2,3);
    zmul = 1 +randn(1)*resize1;%.*obsz(3,3);
    newvert = v1.*repmat([xmul,ymul,zmul],size(v1,1), 1);

    %choose spatial location for point
    arxyzi = arxyz(i,:);
    edloc = [-obsz(1,1)*xmul, 0,0];
    arloc = arxyzi;
    varloc = [0,randn(1,2)]*locsig;

    newloc = repmat(round(arloc + edloc+ varloc), size(v1,1),1);% this must be integer!
%     newloc = repmat([177.1,188,215],size(v1,1),1); 
    verts = newvert+newloc; % fvcombovert = verts doesn't work when newloc is added here, newvert*1.2 is fine
    
    %place a shape in space at a random location
%     p = patch('Vertices', eghandle.vertices+newloc, 'Faces', eghandle.faces);
    p = patch('Vertices', verts , 'Faces', f1);
%     p = isosurface()
    hold on

    % indicate the colour of the shape
    p.FaceColor= pcol(i,:);
    p.EdgeColor= pcol(i,:);

    %rotate the shape according to a distribution 
    %rotate along axis
    rotate(p, [1,0,0], rand(1)*360, newloc(1,:))
%     % rotate off axis
    rotate(p, [0,1,0], randn(1)*rotsig, newloc(1,:))
    rotate(p, [0,0,1], randn(1)*rotsig, newloc(1,:))
        % these rotations aren't being combined into the output array!

        pverts = get(p, 'Vertices');
        pfaces = get(p, 'Faces');
%     px = get(p,'XData');
%     py = get(p,'YData');
%     pz = get(p,'ZData');
%     pV = get(p, '')


if i == 1
    fvcombovert = pverts;
    fvcomboface = pfaces;
%     fvcomboface = f1; fvcombovert = verts; 
else
    tfv = fvcombovert;
    tff = fvcomboface;
    ntfv = length(tfv);
    fvcombovert = [tfv; pverts];
    fvcomboface = [tff; pfaces+ntfv];

end
    
end

% neaten up the plot
looks(gcf)

                             if saveon(2)==1
                            
                                % creates name for the file - numbers are values as identified
                                %     num2str(RADIUS*rad); % radius length
                                %     num2str(PITCHNO*zsc); % this is the pitch
                                %     num2str(reso); % resolution, ie. chunks in one radius
                                %     num2str(t(end)-t(1)/pi); % number of pitches
                                %     num2str(diz*einterval*zsc); % this is the thickness
                                %     num2str((PITCHNO*zsc)/(RADIUS*rad)); % this is the aspect ratio
                                    pitchN = PITCHNO;%round((t(end)-t(1))/pi);
                                    twritename = strcat('hel-ar',num2str(numsh),'-r1',num2str(einterval),'-sz',num2str(rad),'-N',num2str(pitchN),'-AR',num2str(ar),'-THK',num2str(thk));
                                    twritename(twritename=='.')='_';
                            
                            % adds time-date stamp in order not to write over subsequent runs
                                    totname = DT4filename;
                                    twritename=strcat(twritename, totname);
                            
                            % saves stl file to homefol
                                    cd(homefol)
                            %         FV1 = triangulation(f1,v1);
%                             trangstl = triangulation(fv_combined.faces,fv_combined.vertices );
                                
                                    stlwrite(strcat(twritename,'.stl'), fvcomboface,fvcombovert);
                                   
%                                     stlwrite(strcat(twritename,'.stl'), f1,v1);                   % this works

                                    %                                     stlwrite(trangstl, 'eg.stl', 'binary');

                             end


end

%function to find size of 3d object
function obsz = obsz3(eghandle)
    for i = 1:3
        [jx,j] = min(eghandle(:,i));
        [kx,k] = max(eghandle(:,i));
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