%%% RUN THIS SCRIPT TO FULLY REPLICATE FIGURE 14
[V,F] = load_mesh('../data/wwi-plane-triangle-smooth.obj');
[V,F] = meshfix_components(V,F);
[~,C] = connected_components(F);

top = [21 19 20 18 17 9];
bottom = setdiff(1:max(C),top);
%
[VV,TT,FF] = tetgen_components(V,F,'Flags','-q2a0.0001');
[cV,cT,cF,cC,n] = combine_meshes(VV,TT,FF);


off = [3.5 0 0];
R = axisangle2matrix([1 0 0],pi/2)*axisangle2matrix([0 0 1],-pi/3.0);
b = [0.03729 -0.8345 -1.473];
a = [0.03 -0.2869 1.435;[-2 -0.6 -1;2 -0.6 -1]];
rb = 0.15;
ra = 0.15;
bc_fun = @(X) nan_if(normrow(X-b)>rb & min(pdist2(X,a),[],2)>ra)+ ...
  (normrow(X-b)<=rb)*1+ (min(pdist2(X,a),[],2)<=ra).*0;
tsurf(cF,cV*R,'CData',bc_fun(cV),'Edgealpha',0.2);axis equal;view(2)

bVV = VV(bottom);
bTT = TT(bottom);
bFF = FF(bottom);

cZZ = overlap_poisson(VV,TT,bc_fun, @(X) zeros(size(X,1),1), ...
  'K',2,'Solver','quadprog','LowerBound',0,'UpperBound',1);
bZZ = overlap_poisson(bVV,bTT,bc_fun, @(X) zeros(size(X,1),1), ...
  'K',2,'Solver','quadprog','LowerBound',0,'UpperBound',1);

% From here on it's just plotting


[cV,cT,cF,cC,n,cZ] = combine_meshes(VV,TT,FF,cZZ);
E = sharp_edges(cV,cF);
[cCF,I] = cut_edges(cF,E);
cCV = cV(I,:);
cCZ = cZ(I,:);
cCC = cC(I);
[cCV,IM,J] = remove_unreferenced(cCV,cCF);
cCF = IM(cCF);
cCZ = cCZ(J,:);
cCC = cCC(J);
cAO = [];

[bV,bT,bF,bC,n,bZ] = combine_meshes(bVV,bTT,bFF,bZZ);
[bCF,I] = cut_edges(bF,E);
bCV = bV(I,:);
bCZ = bZ(I,:);
bCC = bC(I);
[bCV,IM,J] = remove_unreferenced(bCV,bCF);
bCF = IM(bCF);
bCZ = bCZ(J,:);
bCC = bCC(J);
bAO = [];


if isempty(cAO)
  cAO = ambient_occlusion(cCV,cCF,cCV,per_vertex_normals(cCV,cCF),1000);
  csAO = laplacian_smooth(cCV,cCF,'cotan',[],0.01,'implicit',cAO);
end
if isempty(bAO)
  bAO = ambient_occlusion(bCV,bCF,bCV,per_vertex_normals(bCV,bCF),1000);
  bsAO = laplacian_smooth(bCV,bCF,'cotan',[],0.01,'implicit',bAO);
end


clf;
hold on;
t = {};

caxis([0.01 0.99]);
%colormap(isolines_map((flipud(cbrewer('RdYlBu',10)))));
colormap(((flipud(cbrewer('RdYlBu',10)))));
CM = cbrewer('Set1',(max(cCC)));
rng(2);

cM = sparse(repmat(cCC,1,3),repmat(1:3,size(cCV,1),1),cCV);
bM = sparse(repmat(bCC,1,3),repmat(1:3,size(bCV,1),1),bCV);
CI = knnsearch(full(cM),full(bM));

CM = CM(randperm(end),:);
t{end+1} = tsurf(bCF,bCV*R+0*off,'FaceVertexCData',CM(CI(bCC),:),fphong,fsoft,'EdgeColor','none');
apply_ambient_occlusion(t{end},'SoftLighting',false,'AddLights',false,'AO',matrixnormalize(bsAO).^2,'Factor',0.8);
t{end+1} = tsurf(cCF,cCV*R+2*off,'FaceVertexCData',CM(cCC,:),fphong,fsoft,'EdgeColor','none');
apply_ambient_occlusion(t{end},'SoftLighting',false,'AddLights',false,'AO',matrixnormalize(csAO).^2,'Factor',0.8);

t{end+1} = tsurf(bCF,bCV*R+1*off,'CData',bCZ,'EdgeColor','none',fphong,fsoft);
add_isolines(t{end},'LineWidth',3);
apply_ambient_occlusion(t{end},'SoftLighting',false,'AddLights',false,'AO',matrixnormalize(bsAO).^2,'Factor',0.8);
t{end+1} = tsurf(cCF,cCV*R+3*off,'CData',cCZ,'EdgeColor','none',fphong,fsoft);
add_isolines(t{end},'LineWidth',3);
apply_ambient_occlusion(t{end},'SoftLighting',false,'AddLights',false,'AO',matrixnormalize(csAO).^2,'Factor',0.8);
hold off;
axis equal;
view(0,15);
camlight;
set(gcf,'Color','w')
set(gca,'Visible','off','pos',[0 0 1 1]);
add_shadow([],light('Position',[.2 -1.0 1],'Style','infinite'));
camproj('persp');
