%%%%% EXAMPLE FILE
clc;
clear;
hold off;



% Define boundary condition and right-hand side of equation
%                   LAPLACIAN(u) = rhs
bc_fun = @(V) zeros(size(V,1),1);
rhs_fun = @(V) -6.*ones(size(V,1),1);

% Load n intersecting meshes (in this case, two)
% This example is with 2D triangle meshes, but you can load your own
% tetrahedral mesh here (FA and FB should be indeces to tetrahedra instead
% of triangles).
[VA,FA] = annulus(200,1.2,'R',2);
[VB,FB] = annulus(200,1,'R',1.8);
VV = {VA,VB};
FF = {FA,FB};

% Call to our method
ZZ = overlap_poisson(VV,FF,bc_fun,rhs_fun,'Solver','min_quad_with_fixed');


% Combine meshes and output in order to plot
F = [FA;FB+size(VA,1)];
V = [VA;VB];
Z = [ZZ{1};ZZ{2}];

% Plot
tsurf(F,[V Z],fsoft,fphong)

