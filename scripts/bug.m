%%% RUN THIS SCRIPT TO FULLY REPLICATE FIGURE 22
% Loading meshes
[V,F] = load_mesh('../data/bug.obj');
[V,F] = meshfix_components(V,F);
[V,F] = qslim_components(V,F,10000);
[V,F] = meshfix_components(V,F);
C = connected_components(F);
% Plot meshes
tsurf(F,V,'CData',C,fphong,fsoft,'EdgeAlpha',0.2);
axis equal;
camlight;
drawnow;

% Tetrahedralize meshes
[VV,TT,FF] = tetgen_components(V,F,'Flags','-q2a1');
[cV,cT,cF,cC,n] = combine_meshes(VV,TT,FF);

[VV,TT,FF] = tetgen_components(V,F,'Flags','-q10a10');
[cV,cT,cF,cC,n] = combine_meshes(VV,TT,FF);

% Call method to find msbk for different values of time
ts = [50 200 800];
cZ = zeros(size(cat(1,VV{:}),1),numel(ts));
for i = 1:numel(ts)
  b = snap_points([-11.39 16.09,12.43],cat(1,VV{:}));
  ZZ = overlap_msbk(VV,TT,ts(i),b);
  [cV,cT,cF,cC,n,cZ(:,i)] = combine_meshes(VV,TT,FF,ZZ);
  
  tsurf(cF,cV,'CData',cZ(:,i),fphong,fsoft,'EdgeColor','none');axis equal;camlight;
  colormap(isolines_map(flipud(cbrewer('RdYlGn',20))));
  drawnow;
end


% From here onwards it's just plotting

cC = connected_components(cT);
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

if isempty(cAO)
  cAO = ambient_occlusion(cCV,cCF,cCV,per_vertex_normals(cCV,cCF),1000);
  csAO = laplacian_smooth(cCV,cCF,'cotan',[],0.001,'implicit',cAO);
end

clf; 
hold on;
CM = cbrewer('Set1',(max(cC)));
rng(14);
CM = CM(randperm(end),:);
R = axisangle2matrix([0 0 1],-pi/2.5);;
t = {};
t{end+1} = tsurf(cCF,cCV*R-0.2*off,'FaceVertexCData',CM(cCC,:),fphong,fsoft,'EdgeColor','none');
  apply_ambient_occlusion(t{end},'SoftLighting',false,'AddLights',false,'AO',matrixnormalize(csAO).^2,'Factor',0.8);
off  = [35 0 0];
for i = 1:size(cZ,2)
  t{end+1} = tsurf(cCF,cCV*R+i*off,'CData',cCZ(:,i),fphong,fsoft,'EdgeColor','none');
  add_isolines(t{end},'LineWidth',1);
  hold on;
  apply_ambient_occlusion(t{end},'SoftLighting',false,'AddLights',false,'AO',matrixnormalize(csAO).^2,'Factor',0.8);
end
hold off;
colormap(flipud(cbrewer('RdYlGn',20)));
%colormap(isolines_map(flipud(cbrewer('RdYlGn',20))));
axis equal;camlight;
camproj('persp');
set(gca,'Visible','off','pos',[0 0 1 1]);
set(gcf,'Color','w')
add_shadow([],light('Position',[.2 -0.2 1],'Style','infinite'));
view(0,47);
