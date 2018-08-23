function [angles, nbroktris] = meshAnglesFromFaceEdgeLen2(el2)

coss = el2*[-1 1 1; 1 -1 1; 1 1 -1]*0.5 ./ sqrt(el2(:, [2 3 1]).*el2(:, [3 1 2]));
bti = any( sqrt(el2)*[-1 1 1; 1 -1 1; 1 1 -1]<-eps, 2 );
nbroktris = sum(bti);
if nbroktris>0, warning( sprintf('%d triangles are broken (violates triangle inequality)', nbroktris) ); end
    
coss( bti, : ) = 1;

coss = max( min( coss, 1 ), -1 );

angles = acos(coss);