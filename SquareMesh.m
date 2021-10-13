function [Coord,Elem,Nb,Db]=SquareMesh(Meshtype,L,N)
switch Meshtype
    case 'Crisscross'
        [Coord,Elem,Nb,Db]=crisscross(L,N);
    case 'Right'
      [Coord,Elem,Nb,Db]=rightmesh(L,N);  
    case 'Reverse'
        [Coord,Elem,Nb,Db]=reversemesh(L,N);
    case 'Delaunay'
       [Coord,Elem,Nb,Db]=dulani(1,N); 
end
[n2ed,ed2el]=edge(Elem,Coord);
[Elem,Coord,Db]=addmidpoints(Elem,Coord,Db,ed2el);
end