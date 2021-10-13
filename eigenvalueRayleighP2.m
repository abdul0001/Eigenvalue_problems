function ev=eigenvalueRayleighP2(uh,Coord,Elem) 
iu=0;
iflux=0;
% uh1=uh.^2;
[rspts,qwgts]=Gausspoints(4);
for j=1:size(Elem,1)
     nodes=Elem(j,:);
     [AK,MK]=localmassP2(Coord,nodes,rspts,qwgts);
    C=uh(Elem(j,:))*uh(Elem(j,:))';
    MN1=C.*AK;
    N1=C.*MK;
    iflux=iflux+sum(MN1(:));
    iu=iu+sum(N1(:));
end  
% for j=1:size(Elem,1)
%      iu=iu+1/2*det([1 1 1; Coord(Elem(j,:),:)'])*sum(uh1(Elem(j,:)))/3;
% end
ev=iflux/iu;
end