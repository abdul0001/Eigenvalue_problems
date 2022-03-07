function [Coord,Elem,Nb,Db]=LshapeDomain(Meshtype,N)
%  Meshtype='Mixed';
% N=8;
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
        [Coord,Elem,Nb,Db]=Lmixed(N);
    case 'Delaunay'
       [Coord,Elem,Nb,Db]=dulani(2,1,N); %% 1 is side length
end
% save('mesh/LcrissN64.mat','Coord','Elem','Nb','Db');
% TR=triangulation(Elem,Coord);
% triplot(TR)
end