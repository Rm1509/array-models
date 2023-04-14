function [f1,v1] = makerod(ar, thk, homefol, saveon, r, einterval, options )
        
    %% rods: central radius  a  and lateral radius  r. 
%     ar=4;
%     einterval=0.05;
%     r=100;
%     saveon=0;

    length=r*ar;
    reso= 1/einterval; 
    
    height=length/2;
    zh=height;
    zhs=height/r;
    a1 = r-r/2;
    
    theta = pi*(0:2*reso)/reso;
    phi   = 2*pi*(0:reso)'/reso;
    phi   = pi*(-reso/2:reso/2)'/reso;
    
    x = [zeros(size((a1 + r*cos(phi(1))*sin(theta))));(a1 + r*cos(phi))*cos(theta);zeros(size((a1 + r*cos(phi(1))*sin(theta))))];%*0.2
    y = [zeros(size((a1 + r*cos(phi(1))*sin(theta))));(a1 + r*cos(phi))*sin(theta);zeros(size((a1 + r*cos(phi(1))*sin(theta))))];%*0.5
    
    
    % z = [r*sin(phi)*ones(size(theta))*zhs;flipud(r*sin(phi)*ones(size(theta)))*zhs;r*sin(phi(1))*ones(size(theta))*zhs];
    z = [r*sin(phi(1))*ones(size(theta))*zhs;r*sin(phi)*ones(size(theta))*zhs;r*sin(phi(end))*ones(size(theta))*zhs];
    
    z(z>zh)=zh;
    z(z<-zh)=-zh;
    [polt,polr,polz] = cart2pol(x,y,z);
    polr(polr>r)=r;
    % polr(polr<inner)=inner;
    [x,y,z] = pol2cart(polt, polr,polz);
    figure
    
                    % %     this is how you add noise, in case you want to
                    %     noises = 0.2;
                    %     xnoise1 = rand(size(x,1)-1, size(x,2))*noises;
                    %     xnoise2 = [xnoise1;-xnoise1(1,:)];
                    %     ynoise1 = rand(size(x,1)-1, size(x,2))*noises;
                    %     ynoise2 = [ynoise1;-ynoise1(1,:)];
                    %     znoise1 = rand(size(x,1)-1, size(x,2))*noises;
                    %     znoise2 = [znoise1;-znoise1(1,:)];
                    %     % h = surf(z,y,x,'FaceAlpha',0.5);
                    %     h = surf(z+xnoise2,y+ynoise2,x+znoise2);
    
%     h=surf(x,y,z);
    h=surf(z,y,x);
    alpha 0.5
        xlabel('x')
        ylabel('y')
        zlabel('z')
    axis equal
    [f1,v1,c] = surf2patch(h,'triangles');
    if saveon==1
    twritename = strcat('roundring-outer',num2str(r),'-inner',num2str(inner),'-res',num2str(reso),'-AR',num2str(ar),'curve',num2str(r));
    twritename(twritename=='.')='_';
    cd(homefol)
        stlwrite(strcat(twritename,'.stl'), f1,v1);
    end
end