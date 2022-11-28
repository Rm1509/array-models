function handle1 = makeHelicoid(AR, PITCHNO, thk, einterval, rad, homefol, saveon)

%                   makeHelicoid(3,  1,      [2,2,3], 0.05,   0.05, homefol, 1)
% homefol = 'C:\Users\vg19882\OneDrive - University of Bristol\Documents\bloom\cerinthe Model\matlab workings for shape';
% homefol = 'C:\Users\Rox\OneDrive - University of Bristol\Documents\bloom\cerinthe Model\matlab workings for shape';
    
%     saveon=0;
    
    ar=AR;%ar=3;%for ar=[1,2,3,4,5,6,7];% aspect ratio
    PITCH=PITCHNO; %PITCH = 1; % number of pitches

%     thk = [2,2,3]; % this is the thickness to bulk out the material between points

%     rad = 0.05; % this is the scalar on the radius (R is 1)
%     einterval = 0.05; % resolution
    R = 1.2; % this is the size of the meshgrid
    RADIUS = 1; % this is the unit (unscaled size for the radius)


            xsc = rad;
            ysc = rad;
            zsc = ar*rad; % zsc is the pitch
            diz = ceil(thk(3)/ar); % diz x interval is the thickness
            diy=thk(2);
            dix = thk(1);
    
% this was developed from the notes made in plottingHelicoid.m

% cd(homefol)
% load('pxyz__helicoidal_points_unit_res49.mat')
% load('pxyz__helicoidal_points_unit_SinglePitch.mat')
    %% 2. making a helicoid point space THIS WORKS! TAKES A little while- load below!!
% make xyz points
                            % these are set above
                            % einterval = 0.05;
                            % R = 1.2;
[px,py,pz] = meshgrid(-R:einterval:R);
pxyz = [px(:),py(:),pz(:)];
                            % these are set above
                            % RADIUS = 1;
                            % PITCH = 1;
reso = RADIUS/einterval;

% make xyz function points
t = -PITCH*pi:einterval:PITCH*pi;
t = -PITCH*pi/2:einterval:PITCH*pi/2; % single pitch
u= -RADIUS:einterval:RADIUS;

v = nan(size(px));

hx = u'*cos(t);
hy = u'*sin(t);
hz = (1/3.142)*t ;
hxyz1 = cat(3, hx, hy,ones(size(hx,1),1)*hz );

isbd1 = false(size(hxyz1,1),size(hxyz1,2));
isbd1([1 end], :)=true;
isbd1(:, [1 end])=true;

hxyz = reshape(hxyz1, [],3);
n = size (hxyz1,1);
%% 2a this is the bit from above that takes ages
%%% THIS IS THE FUNCTION THAT TAKES AGES
    k = dsearchn(pxyz,hxyz);
    % k is the valuable thing here - do not delete
    
%     pxyz3d = reshape(pxyz, 161,161,161,3);

%  xxxxxxxxx HERE'S A FIGURE xxxxxxxxxxxxxxxxxx
figure; plot3(pxyz(k,1), pxyz(k,2),pxyz(k,3),'o')
%  xxxxxxxxxxxxxxx END FIGURE xxxxxxxxxxxxxxxxxx
%% 3. here you could optionally LOAD THE THING THAT TAKES AGES ABOVE, rather than remaking it
% NOTHING HERE
%% 4:6 POST 3 (OR 2) this makes a structure from surface, thickens, stretches and SAVES with labels
%---------------------------------------------------------------------------
% 4. this plots the logical structure from the point cloud as a surface
%---------------------------------------------------------------------------
% for i = 1: n
    ttf = false(size(px));
    ttf(k) = true;
                                % THESE ARE SET AT THE BEGINNING; BUT THEY COULD BE SET HERE
                                %     ar=3;%for ar=[1,2,3,4,5,6,7];% aspect ratio
                                %     rad = 0.05; % this is the scalar on the radius (R is 1)
                                %     xsc = rad;
                                %     ysc = rad;
                                %     zsc = ar*0.05; % zsc is the pitch
    dend = [0,0,0];
%     [f1,v1] = isosurface([0:xsc/(size(px,1)-dend(1)):xsc]',[0:ysc/(size(px,2)-dend(2)):ysc]',[0:zsc/(size(px,3)-dend(3)):zsc]', ttf(1:end-dend(1),1:end-dend(2),1:end-dend(3)));
    [f1,v1] = isosurface(squeeze(px(1,1:end-dend(1),1))*xsc, squeeze(py(1:end-dend(2)*ysc,1,1)), squeeze(pz(1,1,1:end-dend(3)))*zsc,ttf);

%  xxxxxxxxx HERE'S A FIGURE xxxxxxxxxxxxxxxxxx
    figure
isosurface(squeeze(px(1,1:end-dend(1),1))*xsc, squeeze(py(1:end-dend(2),1,1))*ysc, squeeze(pz(1,1,1:end-dend(3)))*zsc,ttf)
%  xxxxxxxxxxxx END FIGURE xxxxxxxxxxxxxxxxxx
'4'
%---------------------------------------------------------------------------
% 5. this thickens the surface points into a defined volume
%---------------------------------------------------------------------------
% these are set above
% diz = ceil(3/ar); % diz x interval is the thickness
%     diy=2;
%     dix = 2;
%     

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
    %
    % this was a surface plot of the shape
    
%     isosurface([0:xsc/160:xsc]',[0:ysc/160:ysc]',[0:zsc/158:zsc]', ttf2(:,:,1:end-2));
% ttf2(:,:,159:end) = false;

   [f1,v1] = isosurface(squeeze(px(1,1:end-dend(1),1))'*xsc, squeeze(py(1:end-dend(2),1,1))*ysc, squeeze(pz(1,1,1:end-dend(3)))*zsc,ttf2);
%  xxxxxxxxx HERE'S A FIGURE xxxxxxxxxxxxxxxxxx
    figure
             isosurface(squeeze(px(1,1:end-dend(1),1))*xsc, squeeze(py(1:end-dend(2),1,1))*ysc, squeeze(pz(1,1,1:end-dend(3)))*zsc,ttf2)
%  xxxxxxxxxxxxxxxx END FIGURE xxxxxxxxxxxxxxxxxx
'5'
%---------------------------------------------------------------------------
% 6. turns shape to x axis and saves it with labels
%---------------------------------------------------------------------------
 
num2str(RADIUS*rad); % radius length
num2str(PITCH*zsc); % this is the pitch

num2str(reso); % resolution, ie. chunks in one radius
num2str(t(end)-t(1)/pi); % number of pitches
num2str(diz*einterval*zsc); % this is the thickness
num2str((PITCH*zsc)/(RADIUS*rad)); % this is the aspect ratio

% ttf3 = permute(ttf2(:,:,size(ttf2,2)/2:end),[1,3,2]);
% ttf3 = permute(ttf2,[1,3,2]);

fvcoords = [squeeze(px(1,1:end-dend(1),1))'*xsc,squeeze(py(1:end-dend(2),1,1))*ysc,squeeze(pz(1,1,1:end-dend(3)))*zsc];
fvcoords2 = [fvcoords(:,3),fvcoords(:,1), fvcoords(:,2)];

%  xxxxxxxxx HERE'S A FIGURE xxxxxxxxxxxxxxxxxx
figure; 
handle1 = isosurface(fvcoords2(:,1), fvcoords2(:,2), fvcoords2(:,3), permute(ttf2,[1,3,2]));
[f1,v1] = isosurface(fvcoords2(:,1), fvcoords2(:,2), fvcoords2(:,3), permute(ttf2,[1,3,2]));

% stlwrite(strcat('isosurface_helicoid1_',num2str(diz),'_xsc','.stl'), f1,v1);

    %     stlwrite(FILE, X, Y, Z)
    xlabel('x')
    ylabel('y')
    zlabel('z')

% xlim([-xsc xsc]*2.1/2)
% ylim([-ysc ysc]*2.1/2)
% zlim([-zsc zsc]*2.1/2)
% end
 %  xxxxxxxxxxxxxx END FIGURE xxxxxxxxxxxxxxxxxx

pitchN = round((t(end)-t(1))/pi);

%% this just saves it
if saveon==1

    % this is a name that saves information about the run
twritename = strcat('helicoid1-reso',num2str(reso),'-N',num2str(pitchN),'-AR',num2str(ar),'-THK',num2str(diz*zsc/reso));
twritename(twritename=='.')='_';

    % this adds the time and date in order not to write over subsequent runs
    dtdt=datetime;
 dty =yyyymmdd(dtdt); dty = num2str(dty);dty= dty(3:end);
[h,m,s] = hms(dtdt);
            hh = sprintf('%02d',h);
            mm = sprintf('%02d',m);
            ss = sprintf('%02d',round(s));
            totname = strcat('_DT',dty,'-',hh,mm,ss);
twritename=strcat(twritename, totname);

% this saves
cd(homefol)
FV1 = triangulation(f1,v1);
% stlwrite(strcat(twritename,'.stl'), FV1);

% stlwrite(FV1, strcat(twritename,'.stl'));
    stlwrite(strcat(twritename,'.stl'), f1,v1);
end
%%

    '6'
end
    