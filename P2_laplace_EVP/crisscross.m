function [Coord,Elem,Nb,Db]=crisscross(L,n)
h=L/n;
Elem=zeros(4*n*n,3);
Coord=zeros(2*n*n+2*n+1,2);
Db=zeros(4*n,2);
Nb=[];
N=(n+1)*(n+1);
for i=1:n+1
    for j=1:n+1
        Coord((i-1)*(n+1)+j,:)=[(j-1)*h (i-1)*h];
    end 
end
h1=h/2;
for i=1:n
    for j=1:n
        N=N+1;
         Coord(N,:)=[(2*j-1)*h1 (2*i-1)*h1];
         k=4*(i-1)*n+4*(j-1);
         l=(n+1)*(i-1);
        Elem(k+1,:)=[N l+j l+j+1];
        Elem(k+2,:)=[N l+j+1 l+j+n+2];
        Elem(k+3,:)=[N l+j+n+2 l+j+n+1];
        Elem(k+4,:)=[N l+j+n+1 l+j];
    end
end
for j=1:n
    Db(j,:)=[j j+1];
    Db(j+n,:)=[j*(n+1) j*(n+1)+n+1];
    Db(j+2*n,:)=[(n+1)*(n+1)-(j-1) (n+1)*(n+1)-j];
    Db(4*n-j+1,:)=[(n+1)*j+1 (n+1)*(j-1)+1];
end
TR=triangulation(Elem,Coord);
triplot(TR)
end