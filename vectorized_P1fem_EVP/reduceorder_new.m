function[uh,lambda,np] =reduceorder_new(A,M,S,Coord,Db,ERR) 
% A-stiffness matrix, M-Mass matrix, S-Snapshot matrix, ERR- Tolarence 
% Coord- matrix containd coordinates of nodes, Db- Dirichlet bdy
C=S'*S;
ns=size(S,2);
mk=min(ns,ns);
[U,D]=eigs(C,mk);
d=diag(D);
err1=0;
j=0;
while err1<1-ERR
    j=j+1;
    err1=sum(d(1:j))/sum(d(:));
end
np=j+1;
V=(S*U(:,1:np))./sqrt(d(1:np))';
lambda=1;
rdt=1/0.1;
FullNodes=[1:size(Coord,1)];
FreeNodes=setdiff(FullNodes,unique(Db));
phik=zeros(length(FullNodes),1);
phik(FreeNodes(16:end),1)=1;
psi=V'*phik;
phi0=phik;
A0=A;
A=A+rdt*M;
err=1;
 while err>1e-8
    M1=M*(rdt+lambda);
    bk=V'*M1*V*psi;
    psi=(V'*A*V)\bk;
    phik=V*psi;
    err=norm(phik-phi0)/norm(phik);
    phi0=phik;
    lambda=(phik'*A0*phik)/(phik'*M*phik);
 end
uh=phik;
end