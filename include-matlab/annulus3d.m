function [VV,TT]=annulus3d(n,RH,RL)
    [VH,FH] = lloyd_sphere(n);
    [VL,FL] = lloyd_sphere(n*RH/RL);
    VV = [RH.*VH;RL.*VL];
    FF = [FH;FL+size(VH,1)];
    max_area = 2*avgedge(VV,FF)^3/(6*sqrt(2));
    [VV,TT,FF] = tetgen(VV,FF,'Flags',sprintf('-q1.2a%0.17f',max_area));
    barycentres = .25*(VV(TT(:,1),:)+VV(TT(:,2),:)+VV(TT(:,3),:)+VV(TT(:,4),:));
    I = find(normrow(barycentres)<RL);
    TT(I,:)=[];
    [RV,IM,J] = remove_unreferenced(VV,TT);
    TT = IM(TT);
    VV = RV;
end