%%% RUN THIS SCRIPT TO FULLY REPLICATE FIGURE 4
[V,F] = load_mesh('../data/D00558.ply');
[SVf,SFf] = signed_distance_isosurface(V,F,'GridSize',300,'Level',0.0001,'SignedDistanceType','pseudonormal','ContouringMethod','marching_cubes');
[SVc,SFc] = signed_distance_isosurface(V,F,'GridSize',50,'Level',0.01,'SignedDistanceType','pseudonormal','ContouringMethod','marching_cubes');
C = connected_components(F);
C = C(F(:,1));
R = axisangle2matrix([0 1 0],-pi/2)* ...
  axisangle2matrix([1 0 0],0.08)* ...
  axisangle2matrix([0 0 1],-0.5);
t = {};
off = [1 0 0];
clf;
hold on;
t{end+1} = tsurf(F,V*R+0*off,'EdgeAlpha',0,fsoft,'FaceVertexCData',interp1([max(C);1],[blue+0.7*(1-blue);blue],C));
t{end+1} = tsurf(SFc,SVc*R+1*off,'EdgeAlpha',0,fsoft,'FaceVertexCData',repmat(orange,size(SFc,1),1));
red = [0.9 0.2 0.2];
t{end+1} = tsurf(SFf,SVf*R+2*off,'EdgeAlpha',0,fsoft,'FaceVertexCData',repmat(red,size(SFf,1),1));
hold off;

axis equal;
camproj('persp');
set(gcf,'Color','w');
l = light('Position',[0.2 -0.6 1],'Style','infinite');
s = add_shadow([],l,'Fade','infinite','Color',0.9*get(gcf,'Color'),'BackgroundColor',get(gcf,'Color'));
view(0,20);
camlight;
set(gca,'Visible','off','Position',[0 0 1 1]);
