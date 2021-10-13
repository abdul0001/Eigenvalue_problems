function [Elem2,Coord,Db2]=addmidpoints(Elem,Coord,Db,ed2el)
N=size(Coord,1);
Elem2=[Elem zeros(size(Elem))];
Db2=zeros(2*size(Db,1),2);
cnt=0;
for j=1:size(ed2el,1)
   N=N+1;
   Coord(N,:)=1/2*Coord(ed2el(j,1),:)+1/2*Coord(ed2el(j,2),:);
   l=ed2el(j,3);
   k=ed2el(j,4);
   if all(ed2el(j,1:2)==Elem(l,[1 2]))|| all(ed2el(j,1:2)==Elem(l,[2 1]))
       Elem2(l,4)=N;
   elseif all(ed2el(j,1:2)==Elem(l,[2 3]))||all(ed2el(j,1:2)==Elem(l,[3 2]))
       Elem2(l,5)=N;
   else
      Elem2(l,6)=N; 
   end
   if k==0
       cnt=cnt+1;
       Db2(2*cnt-1,:)=[ed2el(j,1) N];
       Db2(2*cnt,:)=[N ed2el(j,2)];
   end
   
   if k~=0
       if all(ed2el(j,1:2)==Elem(k,[1 2]))|| all(ed2el(j,1:2)==Elem(k,[2 1]))
       Elem2(k,4)=N;
   elseif all(ed2el(j,1:2)==Elem(k,[2 3]))||all(ed2el(j,1:2)==Elem(k,[3 2]))
       Elem2(k,5)=N;
   else
      Elem2(k,6)=N; 
   end    
   end
end  %%% for loop end here

end