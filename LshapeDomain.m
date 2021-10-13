function [Coord,Elem,Nb,Db]=LshapeDomain(Meshtype,N)
n=log2(N);
switch Meshtype
    case 'Crisscross'
        [Coord,Elem,Nb,Db]=Lcrisscross(N);
    case 'Right'
      [Coord,Elem,Nb,Db]=LInitialMesh(1);
      for j=1:n
          [n2ed,ed2el]=edge(Elem,Coord);
          [Coord,Elem,Db,Nb]=redrefine(Coord,Elem,n2ed,ed2el,Db,Nb);
      end
    case 'Reverse'
        [Coord,Elem,Nb,Db]=LInitialMesh(2);
      for j=1:n
          [n2ed,ed2el]=edge(Elem,Coord);
          [Coord,Elem,Db,Nb]=redrefine(Coord,Elem,n2ed,ed2el,Db,Nb);
      end
      case 'Mixed'
        [Coord,Elem,Nb,Db]=LInitialMesh(3);
      for j=1:n
          [n2ed,ed2el]=edge(Elem,Coord);
          [Coord,Elem,Db,Nb]=redrefine(Coord,Elem,n2ed,ed2el,Db,Nb);
      end
    case 'Delaunay'
       [Coord,Elem,Nb,Db]=dulani(2,N); 
end
[n2ed,ed2el]=edge(Elem,Coord);
[Elem,Coord,Db]=addmidpoints(Elem,Coord,Db,ed2el);
end