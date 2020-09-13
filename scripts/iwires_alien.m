%%% RUN THIS SCRIPT TO FULLY REPLICATE FIGURE 1
[OV,OF] = load_mesh('../data/Big.ply');
OV = (OV-min(OV))./max(max(OV)-min(OV));
[DV,DF] = remove_degenerate_faces(OV,OF,'Epsilon',1e-6);
[DV,DF] = meshfix_components(DV,DF);
% 1/6 --> 300,000 vertices
[VV,TT,FF,I] = tetgen_components(DV,DF,'MaxVolFactor',1/6);
[V,T,F] = combine_meshes(VV,TT,FF);

%top = [-130      -0.0002        49.48];
mid =  [0.27395      0.40528      0.81041];

% Select ball at top of "mast"
source_fun = @(V) normrow(V-mid)<0.051;
tsurf(F,V,'CData',1*source_fun(V),'EdgeColor','none',fphong,fsoft);axis equal;camlight;

% Timing of constraint building for reporting in paper table
OO = cellfun(@(F) boundary_faces(F),TT,'UniformOutput',false);
fprintf('overlap_constraints: %g\n',timeit(@() overlap_constraints(VV,TT,'BoundaryFacets',OO)))

lambda = 1000; % Heat geodesics parameter

ZZ = overlap_heat_geodesic(VV,TT,source_fun,'Lambda',lambda);


% From here on it's just plotting

[V,T,F,C,~,Z] = combine_meshes(VV,TT,FF,ZZ);
% Do before cut
plane_n = normalizerow([0 -1 0.05]);
plane_p = [0.59 0.348 0.0929];
plane = [plane_n,-plane_n*plane_p'];
[SV,SF,~,BC] = slice_tets(V,T,plane);
SZ = BC*Z;

[F,I] = cut_edges(F,sharp_edges(V,F));
V = V(I,:);
C = C(I,:);
Z = Z(I,:);

% do after cut
[PV,PF] = half_space_intersect(V,F,plane_p,-plane_n,'Cap',false);
PI = knnsearch(V,PV);
PZ = Z(PI);
% but they get glued again...
[PF,I] = cut_edges(PF,sharp_edges(PV,PF));
PV = PV(I,:);
PZ = PZ(I,:);

[EV,EF] = load_mesh('../data/iwires-alien-explode-view.obj');
EV = EV/max(EV(:))*1.2;
EV = EV-mean(EV)+mean(V);
EV(:,3) = EV(:,3)-min(EV(:,3))+min(V(:,3));
[EV,EF] = upsample(EV,EF,'Iterations',2);
EC = connected_components(EF);
[EF,I] = cut_edges(EF,sharp_edges(EV,EF));
EV = EV(I,:);
EC = EC(I);
R = axisangle2matrix([0 0 1],0.27);
off = [[0.8 0 0].*[0;1;2;3]];
CMS = cbrewer('Set1',max(C));
rng(0);
CMS = CMS(randperm(end),:);
CM = cbrewer('RdYlBu',23);
CMS(end,:) = CM(1,:);


clf;
tsurf(F,V*R+2*off(1,:),'CData',Z,fphong,'EdgeColor','none',fsoft);
tic;AO = apply_ambient_occlusion([],'SoftLighting',false,'AddLights',false);toc
sAO = laplacian_smooth(V,F,'cotan',[],0.001,'implicit',AO);
tsurf(EF,EV*R+off(1,:),'FaceVertexCData',interp1(1:max(EC),CMS,EC(EF(:,1))),'EdgeColor','none',fsoft);
tic;EAO = apply_ambient_occlusion([],'SoftLighting',false,'AddLights',false);toc
tsurf(PF,PV*R+off(1,:),'CData',PZ,'EdgeColor','none',fsoft,fphong);
tic;PAO = apply_ambient_occlusion([],'SoftLighting',false,'AddLights',false);toc
camlight;
axis equal;
view(0,15); 



%




clf;
t = {};
hold on
ca = [min(Z) max(Z)];
t{end+1} = tsurf(EF,EV*R+off(1,:),'FaceVertexCData',interp1(1:max(EC),CMS,EC(EF(:,1))),'EdgeColor','none',fsoft);
t{end+1} = tsurf(F,  V*R+off(2,:),'FaceVertexCData',interp1(1:max(C),CMS,C),fphong,'EdgeColor','none',fsoft);
t{end+1} = tsurf(F,  V*R+off(3,:),'CData',Z,fphong,'EdgeColor','none',fsoft);
ts       = tsurf(SF,SV*R+off(4,:),'CData',SZ,'EdgeColor','none',fphong,fsoft);
t{end+1} = tsurf(PF,PV*R+off(4,:),'CData',(ca(2)-ca(1))+PZ,'FaceAlpha',1,'Edgealpha',0,fphong,fsoft);
tg       = tsurf(F,V*R+off(4,:),'FaceColor',0.9*[1 1 1],'EdgeColor','none',fsoft,'FaceAlpha',0.15);
%colormap(isolines_map(CM));
colormap(CM);
caxis(ca);
view(0,15); 
camlight;
camproj('persp');
axis equal;
add_isolines(t{3},'LineWidth',2);
add_isolines(ts,'LineWidth',2);

apply_ambient_occlusion(t{1},'SoftLighting',false,'AddLights',false,'AO',EAO,'Factor',0.75);
apply_ambient_occlusion(t{2},'SoftLighting',false,'AddLights',false,'AO',sAO, 'Factor',0.75);
apply_ambient_occlusion(t{3},'SoftLighting',false,'AddLights',false,'AO',sAO, 'Factor',0.75);

colormap([CM;hsv2rgb(rgb2hsv(CM).*[1 0.0 0.9])])
caxis([ca(1),ca(2)+(ca(2)-ca(1))]);
add_isolines(t{4},'LineWidth',1,'Color',0.43*[1 1 1]);

apply_ambient_occlusion(t{4},'SoftLighting',false,'AddLights',false,'AO',PAO,'Factor',0.75);

set(gcf,'Color','w');
set(gca,'Visible','off','pos',[0 0 1 1]);
add_shadow(t,light('Position',[.4 -0.6 1],'Style','infinite'), ...
  'Color',get(gcf,'Color')*0.9, ...
  'BackgroundColor',get(gcf,'Color')*0.98);
hold off;
