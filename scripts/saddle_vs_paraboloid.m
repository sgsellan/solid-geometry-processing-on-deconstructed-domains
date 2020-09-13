%%% RUN THIS SCRIPT TO FULLY REPLICATE FIGURE 5
im = im2double(imread('../data/pringles.png'));
tol = 0.2;
[VA,FA] = bwmesh(im(:,:,1)>0.1,'Tol',tol,'SmoothingIters',20);
[VB,FB] = bwmesh(im(:,:,3)>0.1,'Tol',tol,'SmoothingIters',20);
%% Put in unit box
VB = VB-min(VA);
VA = VA-min(VA);
VB = VB/max(max(VA));
VA = VA/max(max(VA));

for pass = 1:2
  switch pass
  case 1
    exact = @(V) ...
      1.75*(-0.6*(V(:,1)-max(VB(:,1))/2).^2 - 1.5*(V(:,2)-max(VA(:,2)/2)).^2);
  case 2
    exact = @(V) ...
      0.6*(V(:,1)-max(VB(:,1))/2).^2 - 1.5*(V(:,2)-max(VA(:,2)/2)).^2;
  end
  
  [An,II] = overlap_constraints({VA,VB},{FA,FB},'Method','naive');
  [Ad,II] = overlap_constraints({VA,VB},{FA,FB},'Method','dirichlet');
  
  I = find(sum(cell2mat(II),2)==2);
  Al = sparse( ...
    repmat((1:numel(I)-1)',1,2), ...
    [repmat(I(1),numel(I)-1,1) I(2:end)], ...
    repmat([1 -1],numel(I)-1,1), ...
    numel(I)-1,size(VA,1)+size(VB,1));
  
  Ze = [exact(VA);exact(VB)];
  I = speye(size(Ze,1),size(Ze,1));
  Zl = min_quad_with_fixed(I,-2*Ze,[],[],Al,zeros(size(An,1),1));
  Zn = min_quad_with_fixed(I,-2*Ze,[],[],An,zeros(size(An,1),1));
  Zd = min_quad_with_fixed(I,-2*Ze,[],[],Ad,zeros(size(Ad,1),1));
  
  R = axisangle2matrix([0 0 1],-pi/6);
  %R = eye(3);
  clf;
  off = [1.15 0 0];
  alpha = 1.0;
  t = {};
  hold on;
  t{end+1} = tsurf(FA,[VA exact(VA)]*R+0*off,'FaceAlpha',alpha,'EdgeColor','none',fphong,fsoft);
  t{end+1} = tsurf(FB,[VB exact(VB)]*R+0*off,'FaceAlpha',alpha,'EdgeColor','none',fphong,fsoft);
  t{end+1} = tsurf(FA,[VA Zl(1:size(VA,1))             ]*R+1*off,'FaceAlpha',alpha,'EdgeColor','none',fphong,fsoft);
  t{end+1} = tsurf(FB,[VB Zl(size(VA,1)+(1:size(VB,1)))]*R+1*off,'FaceAlpha',alpha,'EdgeColor','none',fphong,fsoft);
  t{end+1} = tsurf(FA,[VA Zn(1:size(VA,1))             ]*R+2*off,'FaceAlpha',alpha,'EdgeColor','none',fphong,fsoft);
  t{end+1} = tsurf(FB,[VB Zn(size(VA,1)+(1:size(VB,1)))]*R+2*off,'FaceAlpha',alpha,'EdgeColor','none',fphong,fsoft);
  t{end+1} = tsurf(FA,[VA Zd(1:size(VA,1))             ]*R+3*off,'FaceAlpha',alpha,'EdgeColor','none',fphong,fsoft);
  t{end+1} = tsurf(FB,[VB Zd(size(VA,1)+(1:size(VB,1)))]*R+3*off,'FaceAlpha',alpha,'EdgeColor','none',fphong,fsoft);
  hold off;
  axis equal;
  view(0,15);
  camproj('persp');
  switch pass
  case 1
    colormap(flipud(isolines_map(cbrewer('YlGnBu',8))));
  case 2
    %colormap((isolines_map(cbrewer('YlOrBr',8))));
    colormap(flipud(isolines_map(cbrewer('YlOrRd',8))));
  end
  set(gca,'Visible','off','pos',[0 0 1 1]);
  %set(gca,'Visible','off');
  set(gcf,'Color','w');
  add_shadow([],light('Position',[2.2 -2 3*(pass-1)+1.5],'Style','local'),'Fade','infinite');
  imwrite(myaa({'raw',4}),sprintf('saddle-vs-paraboloid-%d.png',pass));
end
