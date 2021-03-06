% we are solving u_t-\Delta u=\lambda u in \Omega
%              B.C- u=0 on \Gamma
%
clc;
addpath ../code
dt=.1;
rdt=1/dt;
tol1=1e-7; %%tolarance limit used in FOM
np=10; % number of POD basis
Meshtype='Right'; % 'Crisscross', 'Right', Reverse, Delaunay
for i=2:5
err=1;
N=2^i; % number of subinterval in each side of square, usually 2^k
L=pi; % L is side of the square (o,L)^2

%% Create the mesh of the square (0,L)^2 of tupe 
 [Coord,Elem,Nb,Db]=SquareMesh(Meshtype,L,N); % for square domain (0,L)^2
 
 %% Create Lshape domain (-1,1)^2\[0,1]X[-1,0] 'Crisscross', 'Right',
 %'Reverse', 'Delaunay' and 'Mixed'
% [Coord,Elem,Nb,Db]=LshapeDomain(Meshtype,N); % for square domain

%% Fix the sizes of matrices
K1=zeros(size(Coord,1),50);
FullNodes=[1:size(Coord,1)];
FreeNodes=setdiff(FullNodes, unique(Db));
A=sparse(size(Coord,1),size(Coord,1)); % A is global stiffness matrix
M1=sparse(size(Coord,1),size(Coord,1)); % global load vector

%% Initialize the solution
uk=zeros(length(FullNodes),1);
uk(FreeNodes(1,8),1)=1; %% Initial eigenfunction

lambda=zeros(itmax,1);
lambda_k=1;
lambda0=1;
[rspts,qwgts]=Gausspoints(4); % quadrature rule
%% Assembly of A 
% stima is element stiffness matrices
for j=1:size(Elem,1)
    nodes=Elem(j,:);
    [AK,MK]=localmassP2(Coord,nodes,rspts,qwgts);
    A(nodes,nodes)=A(nodes,nodes)+AK+rdt*MK; 
    M1(nodes,nodes)=M1(nodes,nodes)+MK; 
end
k=0;
%% Iteratively solve the Problem
while err>tol1
    k=k+1;
M=(rdt+lambda_k)*M1;
b=M*uk;
uh=solve(Coord,A,b,Nb,Db);
uk=uh;
%% Saving the snapshoots
if mod(k,4)==0
    K1(:,k/4)=uh;
end
%% Calculate Eigenvalues
  lambda_k=eigenvalueRayleighP2(uh,Coord,Elem);
  err=abs(lambda_k-lambda0);
  lambda0=lambda_k;
%  lambda_k=eigenvalue(uh,Coord,Elem,DbElem);
lambda(k,1)=lambda_k;
end
 [rlambda1,it]=reduceorder(Elem,Coord,Db,A,M1,K1,np);
%% Display the computed solution
figure(2)
show(Coord,Elem(:,1:3),uh)
pause(1)
%% Display Eigenvalues
fprintf('The first eigenvalue of FOM with N=%d subintervals:\n',N);
for j=k:k
    fprintf('%.10f\n',lambda(j,1));
end
fprintf('The first eigenvalue of the corresponding ROM :\n');
for j=it:it
    fprintf('%.10f\n',rlambda1(j,1));
end
fprintf('\n');
end