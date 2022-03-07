function [uh,lambda_k,dof,A,M,K,it,err_vec] = evsolveLaplace(Coord,Elem,Db,Nb) 
nE = size(Elem,1);
 lambda_k=1;
 err_vec=zeros(1,100);
FullNodes=[1:size(Coord,1)];
FreeNodes=setdiff(FullNodes, unique(Db));
dof=length(FullNodes);
uk=zeros(length(FullNodes),1);
uk(FreeNodes(1,15),1)=1;
 %*** First vertex of elements and corresponding edge vectors
 c1 = Coord(Elem(:,1),:);
 d21 = Coord(Elem(:,2),:)-c1;
 d31 = Coord(Elem(:,3),:) -c1;
 %*** Vector of element areas 4*|T|
 area4 = 2*(d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1));
 %*** Assembly of stiffness matrix
 I = reshape(Elem(:,[1 2 3 1 2 3 1 2 3])',9*nE,1);
 J = reshape(Elem(:,[1 1 1 2 2 2 3 3 3])',9*nE,1);
 a = (sum(d21.*d31,2)./area4)';
 b = (sum(d31.*d31,2)./area4)';
 c = (sum(d21.*d21,2)./area4)';
 M=ones(nE,1)*[1/6 1/12 1/12 1/12 1/6 1/12 1/12 1/12 1/6];
 M=(M.*area4)'/4;
 A = [-2*a+b+c;a-b;a-c;a-b;b;-a;a-c;-a;c];
 A = sparse(I,J,A(:));
 M = sparse(I,J,M(:));
 K=[];
err=1;
it=0;
rdt=1/0.1;
 while err>1e-8
it=it+1;
A1=A+rdt*M;
M1=(rdt+lambda_k)*M;
b=M1*uk;
uh=zeros(length(FullNodes),1);
uh(FreeNodes)=A1(FreeNodes,FreeNodes)\b(FreeNodes);
err=norm(uh-uk)/norm(uh);
err_vec(it)=err;
uk=uh;
if mod(it,2)==0
K=[K uh];
end
lambda_k=(uh'*A*uh)/(uh'*M*uh);
 end
end