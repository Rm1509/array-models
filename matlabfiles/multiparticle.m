function [fv_combined, volrec] = multiparticle(homefol,numsh, ar, PITCHNO, thk,  saveon, rad, rotsig, resize1,locsig, xloc, spacer1, layers, shape, grid)
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

% grid='hexag'; % or 'square';% or 

% % % %  % save a small box
% % % % cd('C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\hexagonalArrayNumbers\230224\matfiles')
% % % % load('hel-ar7-r10_05-sz100-N1-AR4-THK2  2  3_DT230224-153113.mat', 'volrec')
% % % % box7 = volrec;  

%how many shapes wanted
if ~exist("numsh", "var")
numsh = 1; % number of particles in square grid
end
if ~exist("rotsig","var")
rotsig =0;%10; % sigma of random rotation around y & z axes independently (deg)
end
if ~exist("resize1", "var")
resize1 = 0;%0.2; % sigma of random particle resize uniform in all 3 dimensions
end
if ~exist("locsig", "var")
locsig = 0;%20; % sigma of particle distribution deviation from grid
end
if ~exist("ar", "var")
ar = 1;%3;
end
if ~exist("PITCHNO" ,"var")
PITCHNO = 1;
end
if ~exist("thk", "var")
thk = [2,2,3];
end
if ~exist("saveon", "var")
saveon = [0,0];%[0,1];
end
if ~exist("rad", "var")
rad = 100;
end
if ~exist("spacer1", "var")
spacer1 = 0;
end
einterval=0.05;
% homefol = 'C:\Users\Rox\OneDrive - University of Bristol\Documents\lumerical beb\multiparticleArrayDev\230112';
allvariables.numsh = numsh;
allvariables.rotsig = rotsig;
allvariables.resize1 = resize1;
allvariables.locsig = locsig;
allvariables.xloc = xloc;
allvariables.ar=ar;
allvariables.PITCHNO = PITCHNO;
allvariables.thk = thk;
allvariables.rad = rad;
allvariables.spacer1 = spacer1;
allvariables.shape = shape;
allvariables.layers = layers;
allvariables.grid=grid;
allvariables.homefol=homefol;

%figure setup
% figure
pcol = parula(numsh);

% this is required if you want to save it,  but is not actually important

%this calls the basic shape, for when the same shape is reproduced multiple
%times
                % #shape
if strcmp(shape,'helicoid') || strcmp(shape,'helicoidp')
[f1,v1] = makeHelicoid(ar, PITCHNO, thk, homefol, saveon(1), rad(1), einterval, shape);
elseif strcmp(shape,'rod')
[f1,v1] = makerod(ar, thk, homefol, saveon(1), rad(1), einterval );
elseif strcmp(shape,'ring')
    if size(rad,2)>1
    rad2 = rad(2);
    else
    rad2=0.6*rad(1);
    end
    [f1,v1]=makering(ar, thk, homefol, saveon(1), rad(1), rad2, einterval );
end
% handle1 = makeHelicoid(3,  1,   [2,2,3], homefol, 0,  100, 0.05);
% ar=3; PITCHNO=1; thk = [2,2,3]; saveon=0;  rad = 100; einterval=0.05; homefol='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumericalbeb\221201 multi tests\matfiles'; 

close(gcf);
%how big is the object
size(v1)
obsz = obsz3(v1);

%make an array
if strcmp(grid,'square')
arxyz = squarray(numsh, obsz, spacer1);
elseif strcmp(grid, 'hex')
arxyz = hexray(numsh, obsz, spacer1);
end

% fv_combined = struct([]);
% fv_combined.vertices = [];
% fv_combined.faces = [];
% to make a multiarray, each one is made in this loop

for lno = 1:layers                      
% layoffset= (lno-1)*2*obsz(1,1) ; % -2 makes the thing go backwards
layoffset= (layers/2 - (lno-1))*2*obsz(1,1);
    for i = 1: numsh % going through each one of the number of shapes
    
    %     resizing
        xmul = 1 +randn(1)*resize1;%.*obsz(1,3);
        ymul = 1 +randn(1)*resize1;%.*obsz(2,3);
        zmul = 1 +randn(1)*resize1;%.*obsz(3,3);
        newvert = v1.*repmat([xmul,ymul,zmul],size(v1,1), 1);
    
        %choose spatial location for point
        arxyzi = arxyz(i,:);
        edloc = [layoffset-obsz(1,1)*xmul , 0,0]; % this makes the centre 0, no matter how many layers
        arloc = arxyzi;
        varloc = [randn(1)*xloc,randn(1,2)*locsig];
    
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
            if lno == 1
                Lfvcombovert = fvcombovert;
                Lfvcomboface = fvcomboface;
            %     fvcomboface = f1; fvcombovert = verts; 
            else
                tfv = Lfvcombovert;
                tff = Lfvcomboface;
                ntfv = length(tfv);
                Lfvcombovert = [tfv; fvcombovert];
                Lfvcomboface = [tff; fvcomboface+ntfv];
            
            end
end % for the layers
        fv_combined.vertices =Lfvcombovert;
        fv_combined.faces = Lfvcomboface;
        fvcombovert=Lfvcombovert;
        fvcomboface=Lfvcomboface;
        % here the size of the object group
        % should be measured and returned
x12= [max(fvcombovert(:,1)), min(fvcombovert(:,1))];
y12= [max(fvcombovert(:,2)), min(fvcombovert(:,2))];
z12= [max(fvcombovert(:,3)), min(fvcombovert(:,3))];
volrec = [ceil(abs(x12(1)-x12(2))), ceil(abs(y12(1)-y12(2))),ceil(abs(z12(1)-z12(2)))];
% volrec = [2*ceil(x12(2)), ceil(y12(1)-y12(2)),ceil(z12(1)-z12(2))]; % this make sa big overestimate

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
                                    twritename = strcat('hel-ar',num2str(numsh),'-r1',num2str(einterval),'-sz',num2str(rad(1)),'-N',num2str(pitchN),'-AR',num2str(ar),'-THK',num2str(thk));
                                    twritename(twritename=='.')='_';
                            
                            % adds time-date stamp in order not to write over subsequent runs
                                    totname = DT4filename;
                                    twritename=strcat(twritename, totname);
                            
                            % saves stl file to homefol
                                    cd(strcat(homefol,'/matfiles'))
                            %         FV1 = triangulation(f1,v1);
                                    trangstl = triangulation(fv_combined.faces,fv_combined.vertices );
                                
                                    stlwrite(trangstl, strcat(twritename,'.stl'));
%                 volrec=box7;      % save a small box
                                   save(strcat(twritename,'.mat'),"volrec","allvariables")
%                                     stlwrite(strcat(twritename,'.stl'), f1,v1);                   % this works

                                    %                                     stlwrite(trangstl, 'eg.stl', 'binary');



                             end


end

%function to make a square array of given size
function arxyz = squarray(numsh, obsz, spacer1)
    if numsh==1
        arxm = -obsz(1,3)/2;arym =0; arzm =0;
    %     arxyz =[obsz(1,3)/2,obsz(2,3)/2,obsz(3,3)/2 ];
    else
        arlen = ceil(sqrt(numsh));
        larlen=arlen-1;
        arym = repmat([-larlen/2:larlen/2]', 1, arlen)*(obsz(2,3)+spacer1); %previously this was 1:arlen, but it makes a grid with 0,0 inthe corner
        arzm = repmat([-larlen/2:larlen/2], arlen, 1)*(obsz(3,3)+spacer1);
    %     arxm = zeros(size(arym))*obsz(1,3)/10; %this makes the base 0, instead of the centre
        arxm = -1*ones(size(arym))*obsz(1,3)/2;%this makes the centre at 0, instead of the base
    end    
    arxyz = [arxm(:), arym(:), arzm(:)];
    figure; scatter3(arxm(:), arym(:), arzm(:))
end

%function to make a hexagonal array of given size
function arxyz = hexray(numsh, obsz, spacer1);
    if numsh==1
            arxm = -obsz(1,3)/2;arym =0; arzm =0;
        %     arxyz =[obsz(1,3)/2,obsz(2,3)/2,obsz(3,3)/2 ];
        else
            arlen = ceil(sqrt((numsh/6)*2));
            base=1*(obsz(2,3)+spacer1);
            arzm=[];
            arym=[];
            for num=1:arlen
                z=zeros(num*6,1);
                y=zeros(num*6,1);
                z(1:6)=base*num*cos(2*pi/6.*(0:5));
                y(1:6)=base*num*sin(2*pi/6.*(0:5));
                if num>1
                    for q=1:num-1
                       start_x=z(2)-q*base;
                       radi0=sqrt(start_x^2+y(2)^2);
                       start_alpha=1/3*pi+pi*1/3*1/(num)*q;
                       z(q*6+1:q*6+6)=radi0*cos(start_alpha+pi/3.*(1:6));
                       y(q*6+1:q*6+6)=radi0*sin(start_alpha+pi/3.*(1:6));
                    end
                end
                arzm=[arzm; z];
                arym=[arym; y];

            end
             arzm=[0;arzm];
             arym=[0;arym];
            arxm = zeros(size(arym));%-1*ones(size(arym))*obsz(1,3)/2;%this makes the centre at 0, instead of the base

        end    
        arxyz = [arxm(:), arym(:), arzm(:)];
        figure; scatter3(arxm(:), arym(:), arzm(:))
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