function [Coord,Elem,Nb,Db]=rightmesh(L,n)
%n=4;
h=L/n;
Elem=zeros(2*n*n,3);
Coord=zeros((n+1)*(n+1),2);
Db=zeros(4*n,2);
Nb=[];
for i=1:n+1
    for j=1:n+1
        Coord((i-1)*(n+1)+j,:)=[(j-1)*h (i-1)*h];
    end 
end

for i=1:n
    for j=1:n
         k=2*(i-1)*n+2*(j-1);
         l=(n+1)*(i-1);
        Elem(k+1,:)=[l+j+n+1 l+j l+j+1];
        Elem(k+2,:)=[l+j+1 l+j+n+2 l+j+n+1];
    end
end
for j=1:n
    Db(j,:)=[j j+1];
    Db(j+n,:)=[j*(n+1) j*(n+1)+n+1];
    Db(j+2*n,:)=[(n+1)*(n+1)-(j-1) (n+1)*(n+1)-j];
    Db(4*n-j+1,:)=[(n+1)*j+1 (n+1)*(j-1)+1];
end
% TR=triangulation(Elem,Coord);
% triplot(TR)
end