%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEM for u_t - div . (u_x,u_y)  = \lambda u   in  \Omega
%%%%             (u_x, u_y). n = u_N  on  \partial\Omega_N
%                               u = 0    on  \partial\Omega_D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc;
clear all;
L=pi; % domain is (0,L)^2
Meshtype='Crisscross'; % Right, Reverse, Delaunay
dof0=10;
lam0=1;
rlam0=1;
for i=3:8 
     N=2^i; % number of sub intervals on each sides of square
  [Coord,Elem,Nb,Db]=SquareMesh(Meshtype,L,N);
  % [Coord,Elem,Nb,Db]=LshapeDomain(Meshtype,N); % for L shape domain
tic;
[uh,lam,dof1,A,M,K,it,err_vec]=evsolveLaplace(Coord,Elem,Db,Nb);
toc;
%ex=9.6397238440219; % exact Eigenvalue for L shape domain
ex=2; %  exact Eigenvalue for Square domain
rate=(log(abs(lam0-ex))-log(abs(lam-ex)))/log(sqrt(dof1/dof0)); % FOM rate of Conv.
lam0=lam;
uh=uh/sqrt(uh'*M*uh); % normalizing the eigenvector
ue=uexact(Coord,Elem);
E=ue-uh;
ERR=E'*M*E; % error between the exact and approximate eigenvector
%ERR=1e-4/10^i;
 K1=K(:,2:2:end); % snapshot matrox 
tic;
[uh1,rlam1,np]=reduceorder_new(A,M,K1,Coord,Db,ERR); % ROM soln
toc;
rrate=(log(abs(rlam0-ex))-log(abs(rlam1-ex)))/log(sqrt(dof1/dof0)); % rate ROM
rlam0=rlam1;
dof0=dof1;
% print no of sub intervals, FOM Eigen soln, ROM eigen soln, CPU time FOM
% CPU time FOM, DOF, no of reduced basis, no of itertion in FOM
fprintf('%d\t%1.10f\t%1.10f\t%1.2f\t%1.2f\t%d\t%d\t%d\n',2^i,lam,rlam1,rate,rrate,dof1,np,it);
end
% plot the error
x=1:1:it;
semilogy(x,err_vec(1:it),'b');
xlabel('Iterations') 
ylabel('Relative error') 
title('Error Vs iterations')
