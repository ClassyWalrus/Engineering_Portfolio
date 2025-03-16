function PD=blob()

PD.DistFunc = @Blob;
PD.InitEdgeLen = 0.1;
PD.BBox = [-15,-15;15,15];
PD.RHS = -2*70e9*(5*pi/180);

PD=PD_torsion(PD,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fd = Blob(p)
[TH,R] = cart2pol(p(:,1), p(:,2));

rth = 4 - cos(TH).*sin(7.*TH)
      
fd = R - rth;  

return