%%% RUN THIS SCRIPT TO FULLY REPLICATE FIGURE 20
[V,F] = load_mesh('../data/pistol.obj');
[V,F] = meshfix_components(V,F);
[VV,TT,FF] = tetgen_components(V,F);
[Aeq,II,outer_bb,inner_bb] = overlap_constraints(VV,TT);

fprintf('#T: %d\n',size(cat(1,TT{:}),1));
fprintf('K: %d\n',numel(TT));
OO = cellfun(@(F) boundary_faces(F),TT,'UniformOutput',false);
fprintf('overlap_constraints: %g\n',timeit(@() overlap_constraints(VV,TT,'BoundaryFacets',OO)))

[aAeq,II,outer_bb,inner_bb] = overlap_constraints(VV,TT, ...
  'Method','naive','Sparsification','none');


[LL,MM] = overlap_operators(VV,TT,II);
n = cumsum([0 cellfun(@(V) size(V,1),{VV{:}})]);
[cV,cT,cF,cC,n] = combine_meshes(VV,TT,FF);

Beq = zeros(size(Aeq,1),1);
aBeq = zeros(size(aAeq,1),1);
dt = 0.003;
zdata = struct('force_Aeq_li',true);
ydata = [];struct('force_Aeq_li',true);
a = 1e1;
b = 1e-5;
L = blkdiag(LL{:});
M = blkdiag(MM{:});
Q = ((1+b*dt)*M + a*dt^2*L);

rhs_fun = @(X) 1*(X(:,1)>0.94) + ...
               -0.1*(X(:,1)<=0.94 & X(:,1)>0.90);
ZZ = cellfun(@(V) rhs_fun(V),VV,'UniformOutput',false);
YY = ZZ;
cZ = cat(1,ZZ{:});
cY = cat(1,YY{:});
cZ0 = cZ;
cY0 = cY;

cZZ = [];
for iter = 1:250
  cZZ(:,iter) = cZ;
  if mod(iter,5) == 0
    R = axisangle2matrix([1 0 0],pi);
    clf;
    off = [2 0 0];
    hold on;
    tsurf(cF,cV*R+1*off,'CData',cZ,'EdgeColor','none',fsoft,fphong)
    hold off;
    %caxis(max(abs(cZ))*[-1 1]);
    caxis(0.3*[-1 1]);
    axis equal;
    view(0,15);
    camlight;
    colormap(isolines_map(flipud(cbrewer('RdYlBu',11))));
    %colormap((flipud(cbrewer('RdYlBu',256))));
    colorbar
    drawnow;
    %figgif('pistol.gif');
  end

  B = -M*((2-b*dt)*cZ - cZ0 + dt*0);
  cZ0 = cZ;
  tic;
  [cZ,zdata] = min_quad_with_fixed(Q,2*B,[],[],Aeq,Beq,zdata);
  fprintf('min_quad_with_fixed: %g\n',toc);
  
end

[cV,cT,cF,cC,n] = combine_meshes(VV,TT,FF);
cC = connected_components(cT);
E = sharp_edges(cV,cF);
[cCF,I] = cut_edges(cF,E);
cCV = cV(I,:);
cCZZ = cZZ(I,:);
cCC = cC(I);
[cCV,IM,J] = remove_unreferenced(cCV,cCF);
cCF = IM(cCF);
cCZZ = cCZZ(J,:);
cCC = cCC(J);
cAO = [];
if isempty(cAO)
  cAO = ambient_occlusion(cCV,cCF,cCV,per_vertex_normals(cCV,cCF),1000);
  csAO = laplacian_smooth(cCV,cCF,'cotan',[],0.001,'implicit',cAO);
end
colormap((flipud(cbrewer('RdYlBu',256))));
caxis(0.3*[-1 1]);

RR = axisangle2matrix([0 1 0],0.3)*axisangle2matrix([0 0 1],pi*0.28);
t = {};
ss = [10:20:100];
off = [1.0 0 0];
clf;
hold on;
CM = cbrewer('Set1',(max(cC)));
t{end+1} = tsurf(cCF,cCV*R*RR-0.0*off,'FaceVertexCData',CM(cCC,:),fphong,fsoft,'EdgeColor','none');
for si = 1:numel(ss);
  s = ss(si);
  t{end+1} = tsurf(cCF,cCV*R*RR+si*off,'CData',cCZZ(:,s),'EdgeColor','none',fphong,fsoft);
end
caxis(0.3*[-1 1]);
for si = 1:numel(t);
  apply_ambient_occlusion(t{si},'SoftLighting',false,'AddLights',false,'AO',matrixnormalize(csAO).^2,'Factor',0.8);
end
hold off;

axis equal;
camlight;
camproj('persp');
set(gca,'Visible','off','pos',[0 0 1 1]);
set(gcf,'Color','w')
add_shadow([],light('Position',[-.1 -0.3 1],'Style','infinite'));
view(0,15);
