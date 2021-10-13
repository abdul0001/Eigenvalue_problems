function [rlambda,k]=reduceorder(Elem,Coord,Db,A,M1,K,np)
ns=size(K,2); % number of snapshoots
rdt=1/0.01;
itmax=50;
rlambda=zeros(itmax,1);
mean=1/ns*sum(K');
FullNodes=[1:size(Coord,1)];
FreeNodes=setdiff(FullNodes, unique(Db));
phik=ones(length(FullNodes),1);
% phik(FreeNodes(100,130),2)=1;
N=sparse(size(Coord,1),1); % integral value of Basis functions
for i=1:ns
    K(:,i)=K(:,i)-mean';
end
[U,S,V] = svds(K,np);
psi=U\(phik-mean');
% size(psi)
  lambda_k=eigenvalueRayleighP2(phik,Coord,Elem);
  lambda0=lambda_k;
  err=1;
  k=0;
%  lambda_k=reducedevalue(Coord,psi,U,np,mean',N,DbElem);
while err>1e-6
    k=k+1;
    M=M1*(rdt+lambda_k);
    bk=U'*(M-A)*mean'+U'*M*U*psi;
    psi=(U'*A*U)\bk;
    phik=mean'+U*psi;
% lambda_k=eigenvalueRayleigh(phik,Coord,Elem);
%     lambda_k=reducedevalue(Coord,psi,U,np,mean',N,DbElem);
  lambda_k=eigenvalueRayleighP2(phik,Coord,Elem);
  err=abs(lambda0-lambda_k);
  lambda0=lambda_k;
    rlambda(k,1)=lambda_k;
end

end