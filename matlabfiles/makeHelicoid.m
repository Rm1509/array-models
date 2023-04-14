function [f1,v1]= makeHelicoid(ar, PITCHNO, thk, homefol, saveon, rad, einterval,shape )
%                   makeHelicoid(3,  1,   [2,2,3], homefol, 1,  100, 0.05)
% ar=3; PITCHNO=1; thk = [2,2,3]; saveon=0;  rad = 100; einterval=0.05; homefol='C:\Users\Rox\OneDrive - University of Bristol\Documents\lumericalbeb\221201 multi tests\matfiles'; 

% this function makes a helicoid / helix shape, plots graphs and then saves
% an stl object to the homefol. It was developed from the notes made in plottingHelicoid.m

    % ar        is aspect ratio of helix
    % PITCHNO   is the number of repeated half-wavelength pitches
    % thk       is the widths of the point cloud points (in each dimension) needed to ensure the shape
    %                   doesnt have holes
    % homefol   is the location in which the file should be saved
    % saveon    is whether the file should be saved (1) or not (0)
    % rad       is the radius size of the shape (OPTIONAL)
    % einterval is the resolution of the object mesh grid, or 'resolution' % (OPTIONAL)

% -----------------------------------------------------------------------------------------
arguments
    ar double
    PITCHNO double
    thk double
    homefol string
    saveon int32
    rad double
    einterval double
    shape string
%     options.figadd handle
end


% 1. hard coded and default values
    R = 1.2; % this is the size of the meshgrid
    RADIUS = 1; % this is the unit (unscaled size for the radius)
    dend = [0,0,0]; % defining the end of the shape

    if nargin==6; einterval=0.05; end
    if nargin==5; einterval=0.05; rad=100;end

% 2. setting the dimension scales for the object
    xsc = rad;
    ysc = rad;
    zsc = ar*rad; % zsc is the pitch
    diz = ceil(thk(3)/ar); % diz x interval is the thickness 
    diy=thk(2);
    dix = thk(1);
    
% 3. make xyz points
    [px,py,pz] = meshgrid(-R:einterval:R);
    xran =-R:einterval:R;
    yran=-R:einterval:R;
%     zlow = ; zhigh = 
    zran = (-PITCHNO*3.142/2):einterval:(PITCHNO*3.142/2);
    [px,py,pz] = meshgrid(xran,yran,zran);
    pxyz = [px(:),py(:),pz(:)];
    reso = RADIUS/einterval;

% 4. make xyz function points for a helix
    hxyz = shapehelicoid(PITCHNO, einterval, RADIUS);

% 5. This function takes a long time. It matches each of the new function points to a point on the grid 
    k = dsearchn(pxyz,hxyz); 

%  figure; plot3(pxyz(k,1), pxyz(k,2),pxyz(k,3),'o')%  6. xxxxxxxxx Figure showing the location of the function in the grid x

% 7. This initialises the 4th input to isosurface, ttf, ie. the V of X,Y,Z
    ttf = false(size(px));
    ttf(k) = true;

%   figure; isosurface(squeeze(px(1,1:end-dend(1),1))*xsc, squeeze(py(1:end-dend(2),1,1))*ysc, squeeze(pz(1,1,1:end-dend(3)))*zsc,ttf) %  Here's a figure showing the showing the location of the function in the grid

% 8. thickens points to 3D AND CHOOSES BETWEEN HELICOID AND PARABOLA
% HELICOID
if strcmp(shape, 'helicoid')
    ttf2 = volumise(ttf, dix, diy, diz);         % this makes a normal helicoid
elseif strcmp(shape,'helicoidp')
    ttf2 = volubolise(ttf, dix, diy, diz, px, py); % this makes a parabolic surface helicoid
end
%   figure; isosurface(squeeze(px(1,1:end-dend(1),1))*xsc, squeeze(py(1:end-dend(2),1,1))*ysc, squeeze(pz(1,1,1:end-dend(3)))*zsc,ttf2) %   Figure showing the isosurface 


% 9. turns shape to x axis and labels it in a figure 
%     fvcoords = [squeeze(px(1,1:end-dend(1),1))'*xsc,squeeze(py(1:end-dend(2),1,1))*ysc,squeeze(pz(1,1,1:end-dend(3)))*zsc];
%     fvcoords2 = [fvcoords(:,3),fvcoords(:,1), fvcoords(:,2)];

    fvcx = squeeze(px(1,1:end-dend(1),1))'*xsc; fvcx2 = fvcx(:);
    fvcy = squeeze(py(1:end-dend(2),1,1))*ysc; fvcy2 = fvcy(:);
    fvcz = squeeze(pz(1,1,1:end-dend(3)))*zsc; fvcz2 = fvcz(:);
    
                    %  xxxxxxxxx Here's the final shape figure xxxxxxxxxxxxxxxxxx
    figure; 
%     isosurface(fvcoords2(:,1), fvcoords2(:,2), fvcoords2(:,3), permute(ttf2,[1,3,2]));
%     handle1 = isosurface(fvcoords2(:,1), fvcoords2(:,2), fvcoords2(:,3), permute(ttf2,[1,3,2])); % attempt to return the surface as an output in handle1
%     [f1,v1] = isosurface(fvcoords2(:,1), fvcoords2(:,2), fvcoords2(:,3), permute(ttf2,[1,3,2]));
    [f1,v1] = isosurface(fvcz2, fvcx2, fvcy2, permute(ttf2,[1,3,2]));

    xlabel('x'); ylabel('y'); zlabel('z');% xlim([-xsc xsc]*2.1/2) % ylim([-ysc ysc]*2.1/2) % zlim([-zsc zsc]*2.1/2)   
                    %  xxxxxxxxxxxxxx END FIGURE xxxxxxxxxxxxxxxxxx

%  10. saving the stl file. 
    if saveon==1

    % creates name for the file - numbers are values as identified
    %     num2str(RADIUS*rad); % radius length
    %     num2str(PITCHNO*zsc); % this is the pitch
    %     num2str(reso); % resolution, ie. chunks in one radius
    %     num2str(t(end)-t(1)/pi); % number of pitches
    %     num2str(diz*einterval*zsc); % this is the thickness
    %     num2str((PITCHNO*zsc)/(RADIUS*rad)); % this is the aspect ratio
        pitchN = PITCHNO;%round((t(end)-t(1))/pi);
        twritename = strcat('helicoid1-reso',num2str(reso),'-SZ',num2str(rad),'-N',num2str(pitchN),'-AR',num2str(ar),'-THK',num2str(diz*zsc/reso));
        twritename(twritename=='.')='_';

% adds time-date stamp in order not to write over subsequent runs
        totname = DT4filename;
        twritename=strcat(twritename, totname);

% saves stl file to homefol
        cd(homefol)
%         FV1 = triangulation(f1,v1);
        stlwrite(strcat(twritename,'.stl'), f1,v1);
    end

end % end of function

function ttf2 = volumise(ttf, dix, diy, diz)
 % 8. this thickens the surface points into a defined volume / defined
 % thickness. It creates ttf2, which is used in the isosurface
 
    ttf2 = false(size(ttf));%ttf;
    for p = dix+1: size(ttf,1)
       for q = diy+1:size(ttf,2)
           for r = diz+1: size(ttf,3)-diz-1 % this is important to cap off the object
               if ttf(p,q,r)==true
                   ttf2(p,q,r)=true;
                   if r < size(ttf,3)-diz
                       ttf2(p-dix:p+dix,q-diy:q+diy,r:r+diz)=true;
                   else
                       ttf2(p-dix:p+dix,q-diy:q+diy,r:end)=true;
                   end
                   if r > diz
                       ttf2(p-dix:p+dix,q-diy:q+diy,r:r-diz)=true ;
                   else
                       ttf2(p-dix:p+dix,q-diy:q+diy,r:0)=true ;
                   end
               end
           end
       end
    end

end

function ttf2 = volubolise(ttf, dix, diy, diz, px, py)
    ttf2 = false(size(ttf));%ttf;
    for p = dix+1: size(ttf,1)
       for q = diy+1:size(ttf,2)
           rmid2(p,q) = round((px(1,p,1)^2 + py(q,1,1)^2)*5);
           
           for r = diz+1: size(ttf,3)-diz-1 % this is important to cap off the object
               if ttf(p,q,r)==true
                   ttf2(p,q,r)=true;
                   % ADD A SCHERK SURFACE HERE
                   if r < size(ttf,3)-(diz+rmid2(p,q))
                       ttf2(p-dix:p+dix,q-diy:q+diy,r:r+diz+rmid2(p,q))=true;
                       ttf2(p-dix:p+dix,q-diy:q+diy,r:r+diz+rmid2(p,q))=true;
                   else
                       ttf2(p-dix:p+dix,q-diy:q+diy,r:end)=true;
                   end
%                    if r > diz
%                        ttf2(p-dix:p+dix,q-diy:q+diy,r:r-diz)=true ;
%                    else
%                        ttf2(p-dix:p+dix,q-diy:q+diy,r:0)=true ;
%                    end
               end
           end
       end
    end

end


function hxyz = ring()

end

function hxyz = rod()

end
function hxyz = blob()

end
function hxyz = slab()

end
function hxyz = rosette()

end
function hxyz = ladder()

end
function hxyz = striation()

end