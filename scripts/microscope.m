%%% RUN THIS SCRIPT TO FULLY REPLICATE FIGURE 2
[V,F] = load_mesh('../data/microscope.obj');
[V,F] = meshfix_components(V,F);
[SV,SF] = mesh_boolean(V,F,[],[],'union');
[SV,SF] = mesh_boolean(SV,SF,[],[],'union');
try
[TV,TT,TF] = tetgen(SV,SF);
assert(false);
catch e
  % yup, tetgen fails
end
% If you have quartet installed, you can uncomment this to validate our
% result
%[QV,QT,QF] = quartet(SV,SF,0.25);
% writeMESH('../data/microscope-quartet-0.25.mesh',QV,QT,QF);
[QV,QT,QF] = readMESH('../data/microscope-quartet-0.25.mesh');

E = sharp_edges(V,F);
[CF,I] = cut_edges(F,E);
CV = V(I,:);
C = connected_components(F);
CC = C(I);
AO = [];
cAO = [];
QQAO = [];


bc_fun = @(X) nan_if(point_mesh_squared_distance(X,V,F(C==1,:))>1e-3)+1;
rhs_fun = @(X) (point_mesh_squared_distance(X,V,F(C==1,:))<1e-3)*-1;
[VV,TT,FF] = tetgen_components(V,F,'Flags','-q1.2a0.5');
lambda = 10000;
ZZ = overlap_poisson(VV,TT,@(X) nan(size(X,1),1),rhs_fun,'Lambda',lambda,'Solver','min_quad_with_fixed-force_Aeq_li');
%YY = overlap_poisson(VV,TT,@(X) nan(size(X,1),1),rhs_fun,'Lambda',lambda,'Method','naive');
[cV,cT,cF,cC,n,cZ] = combine_meshes(VV,TT,FF,ZZ);
%[ ~, ~, ~, ~,~,cY] = combine_meshes(VV,TT,FF,YY);
lambda = 10000;
QL = -cotmatrix(QV,QT);
QM = massmatrix(QV,QT);
QZ = min_quad_with_fixed(lambda*QL+QM,2*QM*rhs_fun(QV),[],[]);


R = axisangle2matrix([1 0 0],-pi/2);
RR = axisangle2matrix([0 0 1],0.0);

clf;
hold on;
t = {};
off = [80 0 0];
CM = cbrewer('Set1',(max(C)));
t{end+1} = tsurf(CF,CV*R*RR','FaceVertexCData',CM(CC,:),'EdgeAlpha',0.1,fphong,fsoft);

if isempty(AO)
  AO = ambient_occlusion(CV,CF,CV,per_vertex_normals(CV,CF),1000);
  sAO = laplacian_smooth(CV,CF,'cotan',[],0.01,'implicit',AO);
end
apply_ambient_occlusion(t{1},'SoftLighting',false,'AddLights',false,'AO',matrixnormalize(sAO).^2,'Factor',0.8);


E = sharp_edges(cV,cF);
[cCF,I] = cut_edges(cF,E);
cCV = cV(I,:);
cCZ = cZ(I,:);
[cCV,IM,J] = remove_unreferenced(cCV,cCF);
cCZ = cCZ(J,:);
t{end+1} = tsurf(cCF,cCV*R+1*off,'CData',matrixnormalize(cCZ),fsoft,fphong,'EdgeColor','none');


E = sharp_edges(QV,QF);
[QQF,I] = cut_edges(QF,E);
QQV = QV(I,:);
QQZ = QZ(I,:);
[QQV,IM,J] = remove_unreferenced(QQV,QQF);
QQZ = QQZ(J,:);
t{end+1} = tsurf(QQF,QQV*R*RR+2*off,'CData',matrixnormalize(QQZ),fsoft,fphong,'EdgeColor','none');

hold off;
%colormap(isolines_map(flipud(cbrewer('RdYlBu',10))));
colormap(flipud(cbrewer('RdYlBu',20)));
add_isolines(t(2:3),'LineWidth',3);

if isempty(cAO)
  cAO = ambient_occlusion(cCV,cCF,cCV,per_vertex_normals(cCV,cCF),1000);
  csAO = laplacian_smooth(cCV,cCF,'cotan',[],0.01,'implicit',cAO);
end
apply_ambient_occlusion(t{2},'SoftLighting',false,'AddLights',false,'AO',matrixnormalize(csAO).^2,'Factor',0.8);

if isempty(QQAO)
  QQAO = ambient_occlusion(QQV,QQF,QQV,per_vertex_normals(QQV,QQF),1000);
  QQsAO = laplacian_smooth(QQV,QQF,'cotan',[],0.01,'implicit',QQAO);
end
apply_ambient_occlusion(t{3},'SoftLighting',false,'AddLights',false,'AO',matrixnormalize(QQsAO).^2,'Factor',0.8);


axis equal;
view(0,8);
camproj('persp');
set(gca,'Visible','off','pos',[0 0 1 1]);
set(gcf,'Color','w')
add_shadow([],light('Position',[.2 -0.2 1],'Style','infinite'));

camlight;
