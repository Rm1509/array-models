
function hxyz = shapehelicoid(PITCHNO, einterval, RADIUS)%, px)
% this function takes the inputs and makes an array of points that lie on
% the equation for a 'helicoid'
    t = -PITCHNO*pi/2:einterval:PITCHNO*pi/2; % single pitch
    u= -RADIUS:einterval:RADIUS;
    
%     v = nan(size(px));
    
    hx = u'*cos(t);
    hy = u'*sin(t);
    hz = (1/3.142)*t ;
    hxyz1 = cat(3, hx, hy,ones(size(hx,1),1)*hz );
    
    isbd1 = false(size(hxyz1,1),size(hxyz1,2));
    isbd1([1 end], :)=true;
    isbd1(:, [1 end])=true;
    
    hxyz = reshape(hxyz1, [],3);
end