function [Coord,Elem,Nb,Db]=SquareMesh(Meshtype,L,N)
%   Meshtype='Crisscross';
%  L=1;N=4;
switch Meshtype
    case 'Crisscross'
        [Coord,Elem,Nb,Db]=crisscross(L,N);
    case 'Right'
      [Coord,Elem,Nb,Db]=rightmesh(L,N);  
    case 'Reverse'
        [Coord,Elem,Nb,Db]=reversemesh(L,N);
    case 'Delaunay'
       [Coord,Elem,Nb,Db]=dulani(1,L,N); 
end
%  save('DelaunayS512.mat','Coord','Elem','Nb','Db');
% TR=triangulation(Elem,Coord);
% triplot(TR)
% saveas(triplot(TR),'test.png');
end