function uex=uexact(Coord,Elem)
x=Coord(:,1);y=Coord(:,2);
ue=sin(x).*sin(y);
uex=2/pi*ue;
% uex=ue/sqrt(ue'*M*ue);
 trisurf(Elem,x,y,uex,'facecolor','interp')
title ( 'Normalized Eigenfunction')

end