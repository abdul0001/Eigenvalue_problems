%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEM for u_t - div . (u_x,u_y)  = \lambda u   in  \Omega
%%%%             (u_x, u_y). n = u_N  on  \partial\Omega_N
%                               u = 0    on  \partial\Omega_D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc;
clear all;
L=pi;
Meshtype='Crisscross';
dof0=10;
lam0=1;
rlam0=1;
for i=3:8
     N=2^i;
  [Coord,Elem,Nb,Db]=SquareMesh(Meshtype,L,N);

tic;
[uh,lam,dof1,A,M,K,it,err_vec]=evsolveLaplace1(Coord,Elem,Db,Nb);
toc;
%ex=9.6397238440219;
ex=2;
rate=(log(abs(lam0-ex))-log(abs(lam-ex)))/log(sqrt(dof1/dof0));
lam0=lam;
uh=uh/sqrt(uh'*M*uh);
ue=uexact(Coord,Elem);
E=ue-uh;
ERR=E'*M*E;
tic;
[uh1,rlam1,np]=reduceorder_new(A,M,K,Coord,Db,ERR);
toc;
rrate=(log(abs(rlam0-ex))-log(abs(rlam1-ex)))/log(sqrt(dof1/dof0));
rlam0=rlam1;
dof0=dof1;
fprintf('%d\t%1.10f\t%1.10f\t%1.2f\t%1.2f\t%d\t%d\t%d\n',2^i,lam,rlam1,rate,rrate,dof1,np,it);
end
x=1:1:it;
semilogy(x,err_vec(1:it),'b');
xlabel('Iterations') 
ylabel('Relative error') 
title('Error Vs iterations')
