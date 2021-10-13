function [AK,MK]=localmassP2(Coord,nodes,rspts,qwgts)
x=Coord(nodes,1); % node x-coordinates
    y=Coord(nodes,2); % y-
    AK=zeros(6,6); % elements stiffness
    MK=zeros(6,6);
    for q=1:length(qwgts) % quadrature loop
        r=rspts(q,1); % quadrature r-coordinate 
        s=rspts(q,2); % s- 
        [S,dSdx,dSdy,detJ]=Isopmap(x,y,r,s,@P2shapes); 
        wxarea=qwgts(q)*detJ/2; % weight times area 
        AK=AK+(dSdx*dSdx'...
           +dSdy*dSdy')*wxarea; % element stiffness
       MK=MK+(S*S')*wxarea; % 
    end
end
