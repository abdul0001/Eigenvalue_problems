function [Coord,Elem,Db,Nb]=redrefine(Coord,Elem,n2ed,ed2el,Db,Nb)

nt=size(Elem,1); 
ne=size(ed2el,1);
marker=sparse(ne,1);
DbElem=zeros(size(Db,1),3);
for i=1:ne
    inode=ed2el(i,1); enode=ed2el(i,2);
    coord1=Coord(inode,:); coord2=Coord(enode,:);
    nCoord=(coord1+coord2)/2;
    marker(i)=size(Coord,1)+1; 
    Coord(size(Coord,1)+1,:)=[nCoord(1) nCoord(2)];
end


for i=1:nt
    ct=Elem(i,:);
    ce=diag(n2ed(ct([2 3 1 ]),ct([3 1 2])));
    m1=marker(ce(1));  m2=marker(ce(2)); m3=marker(ce(3));
    Elem(i,:)=[m1 m2 m3];
    nt1=size(Elem,1); 
    Elem(nt1+1,:)=[ct(1) m3 m2];  
    Elem(nt1+2,:)=[ct(2) m1 m3];  
    Elem(nt1+3,:)=[ct(3) m2 m1];  
end


%%%% Boundary Edges
if (~isempty(Db))
    nb=size(Db,1);
    for i=1:nb
        base=n2ed(Db(i,1),Db(i,2));
        if (marker(base)>0)
            p=[Db(i,1) marker(base)  Db(i,2)];
            DbElem(i,:)=p;
            Db(i,:)=[p(1) p(2)];
            Db(size(Db,1)+1,:)=[p(2) p(3)];
        end
    end
end


if (~isempty(Nb))
    nb=size(Nb,1);
    for i=1:nb,
         base=n2ed(Nb(i,1),Nb(i,2));
        if (marker(base)>0)
            p=[Nb(i,1) marker(base)  Nb(i,2)];
            Nb(i,:)=[p(1) p(2)];
            Nb(size(Nb,1)+1,:)=[p(2) p(3)];
        end
    end
end
TR=triangulation(Elem,Coord);
triplot(TR)
