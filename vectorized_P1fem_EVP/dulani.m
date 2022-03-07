function [Coord,Elem,Nb,Db]=dulani(N,l,n)
%l=L;

if N==1
 h=l/n;
Nb=[];
g=[2 0 l 0 0 1 0;2 l l 0 l 1 0;
   2 l 0 l l 1 0;2 0 0 l 0 1 0];
[p,e,t] = initmesh(g','hmax',h);
Coord=p';
Elem=t(1:3,:)';
Db=e(1:2,:)';
end
if N==2
h=1/n;
Nb=[];
g=[2 0 1 0 0 1 0;2 1 1 0 1 1 0;
   2 1 -1 1 1 1 0;2 -1 -1 1 -1 1 0
   2 -1 0 -1 -1 1 0;2 0 0 -1 0 1 0];
[p,e,t] = initmesh(g','hmax',h);
Coord=p';
Elem=t(1:3,:)';
Db=e(1:2,:)';
end
%pdemesh(p,e,t);
 TR=triangulation(Elem,Coord);
 triplot(TR)
end