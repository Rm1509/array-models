function hxyz = shaperod(length, einterval, RADIUS)
                        
                        scalar = 10;
                        outer = 210./2;  %5*scalar;
                        r = RADIUS; %20*scalar;   %2.9*scalar; % r<outer-inner
                        reso= 1/einterval; 
                        
                        
                        height=length/2;
                        zh=height;
                        zhs=height/r;
                        a1 = outer-r/2;
                        
                        theta = pi*(0:2*reso)/reso;
%                         phi   = 2*pi*(0:reso)'/reso;
                        phi   = pi*(-reso/2:reso/2)'/reso;
                        
                        
                        x = [zeros(size((a1 + r*cos(phi(1))*sin(theta))));(a1 + r*cos(phi))*cos(theta);zeros(size((a1 + r*cos(phi(1))*sin(theta))))];%*0.2
                        y = [zeros(size((a1 + r*cos(phi(1))*sin(theta))));(a1 + r*cos(phi))*sin(theta);zeros(size((a1 + r*cos(phi(1))*sin(theta))))];%*0.5
                        
                        z = [r*sin(phi(1))*ones(size(theta))*zhs;r*sin(phi)*ones(size(theta))*zhs;r*sin(phi(end))*ones(size(theta))*zhs];
                        
                        z(z>zh)=zh;
                        z(z<-zh)=-zh;
                        [polt,polr,polz] = cart2pol(x,y,z);
                        polr(polr>outer)=outer;

                        [x,y,z] = pol2cart(polt, polr,polz);
                        
                        hxyz=[x,y,z];

                        
end
    








