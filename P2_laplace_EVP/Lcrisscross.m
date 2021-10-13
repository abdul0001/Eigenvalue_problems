function [Coord,Elem,Nb,Db]=Lcrisscross(n)
% n=2;
h=1/n;
Elem=zeros(12*n*n,3);
Coord=zeros((3*n+1)*(n+1)+3*n*n,2);
Db=zeros(8*n,2);
Nb=[];
N=3*(n+1)*(n+1)-2*(n+1);
%%% Coordinate of vertices in [-1,0]X[-1,1]
for i=1:2*n+1
    for j=1:n+1
        k=(i-1)*(n+1)+j;
        Coord(k,:)=[-1+(j-1)*h -1+(i-1)*h];
    end 
end
% Coordinate of vertices in [0,1]X[0,1]
for i=1:n+1
    for j=2:n+1
        Coord(k+(i-1)*n+j-1,:)=[(j-1)*h (i-1)*h];
    end
end

h1=h/2;
for i=1:2*n
    for j=1:n
        N=N+1;
         Coord(N,:)=[-1+(2*j-1)*h1 -1+(2*i-1)*h1];
         k=4*(i-1)*n+4*(j-1);
         l=(n+1)*(i-1);
        Elem(k+1,:)=[N l+j l+j+1];
        Elem(k+2,:)=[N l+j+1 l+j+n+2];
        Elem(k+3,:)=[N l+j+n+2 l+j+n+1];
        Elem(k+4,:)=[N l+j+n+1 l+j];
    end
end
K=8*n*n;
L=2*(n+1)^2-(n+1);
for i=1:n
    for j=1:n
         N=N+1;
         Coord(N,:)=[(2*j-1)*h1 (2*i-1)*h1];
         k=4*(i-1)*n+4*(j-1);
         l=L+n*(i-1);
        if j==1
         l1=(n+1)^2+(i-1)*(n+1);
        Elem(K+k+1,:)=[N l1 l+j];
        Elem(K+k+2,:)=[N l+j l+j+n];
        Elem(K+k+3,:)=[N l+j+n l1+n+1];
        Elem(K+k+4,:)=[N l1+n+1 l1];
        else
        Elem(K+k+1,:)=[N l+j-1 l+j];
        Elem(K+k+2,:)=[N l+j l+j+n];
        Elem(K+k+3,:)=[N l+j+n l+j+n-1];
        Elem(K+k+4,:)=[N l+j+n-1 l+j-1];
        end
        
    end
end


for j=1:n
    Db(j,:)=[j j+1];
    Db(j+n,:)=[j*(n+1) j*(n+1)+n+1];
    Db(j+n*3,:)=[L+j*n L+(j+1)*n];
    Db(j+n*5,:)=[L-(j-1) L-j];
    if j>1
        m=3*(n+1)^2-2*(n+1);
        Db(2*n+j,:)=[L+j-1 L+j];
        Db(4*n+j-1,:)=[m-(j-2) m-(j-1)];
    end
end
for j=1:2*n
    m=2*n*(n+1)+1;
    Db(6*n+j,:)=[m-(j-1)*(n+1) m-j*(n+1)];
end
Db(2*n+1,:)=[(n+1)^2 L+1];
Db(5*n,:)=[L+n*n+1 L];
TR=triangulation(Elem,Coord);
triplot(TR)
end