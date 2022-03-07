function [Coord,Elem,Nb,Db]=Lmixed(n)
%n=32;
h=1/n;

Elem=zeros(6*n*n,3);
Coord=zeros((3*n+1)*(n+1),2);
Db=zeros(8*n,2);
Nb=[];
n1=(2*n+1)*(n+1);
%Coord=zeros(n1,2);
idx=1:n1;
idx=idx';
N=3*(n+1)*(n+1)-2*(n+1);
%%% Coordinate of vertices in [-1,0]X[-1,1]
Coord(idx,:)=[-1+h*(mod(idx-1,n+1)) -1+h*floor((idx-1)./(n+1))];
%Coordinate of vertices in [0,1]X[0,1]
idx1=1:(n+1)*n;
idx1=idx1';
Coord(n1+idx1,:)=[h*(mod(idx1-1,n)+1)  h*floor((idx1-1)./n)];


for i=1:2*n
    for j=1:n
         k=2*(i-1)*n+2*(j-1);
       l=(n+1)*(i-1);
       if i<n+1
        Elem(k+1,:)=[l+j l+j+1 l+j+n+2];
        Elem(k+2,:)= [l+j+n+2 l+j+n+1 l+j];
       else 
        Elem(k+1,:)=[l+j+n+1 l+j l+j+1];
        Elem(k+2,:)=[l+j+1 l+j+n+2 l+j+n+1];
       end
    end
end

K=4*n*n;
L=2*(n+1)^2-(n+1);
for i=1:n
    for j=1:n
         k=2*(i-1)*n+2*(j-1);
         l=L+n*(i-1);
        if j==1
         l1=(n+1)^2+(i-1)*(n+1);
        Elem(K+k+1,:)=[l1 l+j+n l1+(n+1)];
        Elem(K+k+2,:)=[l1 l+j l+j+n];
        else
        Elem(K+k+1,:)=[l+j-1 l+j l+j+n];
        Elem(K+k+2,:)=[l+j-1 l+j+n l+j+n-1];
        end
        
    end
end
%Elem(4*n:end,:)
% %%% Boundary Edges
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
% toc;
% TR=triangulation(Elem,Coord);
% triplot(TR)
end