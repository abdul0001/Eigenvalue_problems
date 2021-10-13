function uh=solve(Coord,A,b,Nb,Db)
FullNodes=[1:size(Coord,1)];
FreeNodes=setdiff(FullNodes, unique(Db));
if (~isempty(Nb))
   for j=1:size(Nb,1)
       b(Nb(j,:))= b(Nb(j,:))+norm(Coord(Nb(j,1),:)-Coord(Nb(j,2),:))*...
                                     u_N(Coord(Nb(j,1),:),Coord(Nb(j,2),:))/2;
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uh=zeros(length(FullNodes),1);
% Dirichlet Conditions
if (~isempty(Db))
    Dbnodes=unique(Db);
    for j=1:size(Dbnodes,1)
       uh(Dbnodes(j),1)=0;
    end
end
b=b-A*uh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solving the linear system
uh(FreeNodes)=A(FreeNodes,FreeNodes)\b(FreeNodes);
end