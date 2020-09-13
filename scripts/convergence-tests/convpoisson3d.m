%%%%%%%%%%%%%%-3D POISSON EQUATION CONVERGENCE TEST-%%%%%%%%%%%%%%%%%%%%%%%
%   This script is designed to exhaustively test convergence of the 3D poisson
%equation solver with dirichlet (outer) boundary conditions, as in
%                           lap u = f,
%               u|_{\partial\Omega} = g,
%where f,g are known functions (rhs_fun and bc_fun respectively) and an
%element discretization of \Omega is not known. Instead, we use element
%discretizations of intersectin domains \Omega_i where
%                    \cup\Omega_i = \Omega.
%   We will test two different 'intersection' constraints: 'naive' and
%'dirichlet', which are in essence equality constraints on the intersection
%or the boundary of the intersection (inner boundaries) of the domains,
%respectively. Mathematically, this amounts to for all j,i=1,...,k,
% 'naive':    
%            u_i|_{\Omega_i\cap\Omega_j}=u_j|_{\Omega_i\cap\Omega_j}
% 'dirichlet':    
%  u_i|_{\partial(\Omega_i\cap\Omega_j)}=u_j|_{\partial(\Omega_i\cap\Omega_j)}
%
%   In most example we also compute the numerical 'ground-truth', where we
%solve on a discretization of the whole \Omega, for reference.
%
%
%
%%%%%%%%%%%%%%%%%%%% LAPLACE EQUATION ON 3D ANNULUS %%%%%%%%%%%%%%%%%%%%%%%
%   We consider \Omega to be a 3D annulus of radii 1 and 2, which we can see
%as the intersecting union of \Omega_1 and \Omega_2, annuli of radii 1 and
%1.8 and 1.2 and 2, respectively. In this case, we solve the equation
%                             lap u = 0
% subject to
%                    u|_{R=1} = u|_{R=2} = 1-1/r
%
%%%%%% Finding analytical solution with separation of variables %%%%%%%%%%%
figure(1)
clear
hold off
phi = @(V) 1-1./normrow(V);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Loop to find and plot errors %%%%%%%%%%%%%%%%%%%%%%%%%%%
bc_fun = phi;
rhs_fun = @(V) zeros(size(V,1),1);
for s=1:5
    [VA,FA] = annulus3d(2^(s+3),2,1.2);
    [VB,FB] = annulus3d(2^(s+3),1.8,1);
    [V,F] = annulus3d(2^(s+3),2,1);
    VV = {VA,VB};
    FF = {FA,FB};
    ZZd = overlap_poisson(VV,FF,bc_fun,rhs_fun,'Method','dirichlet');
    ZZn = overlap_poisson(VV,FF,bc_fun,rhs_fun,'Method','naive');
   % ZZgt = overlap_poisson3({V},{F},bc_fun,rhs_fun);
    av(s) = avgedge([VA;VB],[FA;FB+size(VA,1)]);
    %avgt(s) = avgedge(V,F);
    errord(s) = max(abs(phi([VA;VB])-[ZZd{1};ZZd{2}]));
    errorn(s) = max(abs(phi([VA;VB])-[ZZn{1};ZZn{2}]));
    %errorgt(s) = max(abs(phi(V)-ZZgt{1}));
    plot(log(av),log(errord),log(av),log(errorn),'LineWidth',3)
    legend('dirichlet','naive')
    title('Laplace equation convergence')
    xlabel('log h')
    ylabel('log max error')
    axis equal
    drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% POISSON EQUATION ON 3D ANNULUS %%%%%%%%%%%%%%%%%%%%%%%
%   We consider \Omega to be an annulus of radii 1 and 2, which we can see
%as the intersecting union of \Omega_1 and \Omega_2, annuli of radii 1 and
%1.8 and 1.2 and 2, respectively. In this case, we solve the equation
%                             lap u = 6
%   subject to
%                       u|_{R=1} = u|_{R=2} = r^2-1
%
%%%%%% Finding analytical solution with separation of variables %%%%%%%%%%%
figure(2)
clear
hold off
phi = @(V) (normrow(V).^2)-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Loop to find and plot errors %%%%%%%%%%%%%%%%%%%%%%%%%%%
bc_fun = phi;
rhs_fun = @(V) 0*V(:,1)+6; 
for s=1:7
    [VA,FA] = annulus3d(2^(s+3),2,1.2);
    [VB,FB] = annulus3d(2^(s+3),1.8,1);
    [V,F] = annulus3d(2^(s+3),2,1);
    VV = {VA,VB};
    FF = {FA,FB};
    ZZd = overlap_poisson(VV,FF,bc_fun,rhs_fun,'Method','dirichlet');
    ZZn = overlap_poisson(VV,FF,bc_fun,rhs_fun,'Method','naive');
    ZZgt = overlap_poisson({V},{F},bc_fun,rhs_fun);
    av(s) = avgedge([VA;VB],[FA;FB+size(VA,1)]);
    avgt(s) = avgedge(V,F);
    errord(s) = max(abs(phi([VA;VB])-[ZZd{1};ZZd{2}]));
    errorn(s) = max(abs(phi([VA;VB])-[ZZn{1};ZZn{2}]));
    errorgt(s) = max(abs(phi(V)-ZZgt{1}));
    plot(log(av),log(errord),log(av),log(errorn),...
        log(avgt),log(errorgt),'LineWidth',3)
    legend('dirichlet','naive','groundtruth')
    title('Poisson equation convergence')
    xlabel('log h')
    ylabel('log max error')
    axis equal
    drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% NON RADIALLY-SYMMETRIC CONVERGENCE %%%%%%%%%%%%%%%%%%%%%%%%
%   We proceed now in a similar way, only know we pick any function 'phi'
%well defined on an annulus and compute its laplacian to test the
%convergence of our method for non radially-simmetric functions or boundary
%conditions.
figure(3)
clear
hold off
phi = @(V) V(:,1).*V(:,2).*(V(:,3).^2);
bc_fun = phi;
rhs_fun = @(V) 2.*V(:,1).*V(:,2); %!!!!!!!!!!!!!!!!!!!!
for s=1:7
    [VA,FA] = annulus3d(2^(s+3),2,1.2);
    [VB,FB] = annulus3d(2^(s+3),1.8,1);
    [V,F] = annulus3d(2^(s+3),2,1);
    VV = {VA,VB};
    FF = {FA,FB};
    ZZd = overlap_poisson(VV,FF,bc_fun,rhs_fun,'Method','dirichlet');
    ZZn = overlap_poisson(VV,FF,bc_fun,rhs_fun,'Method','naive');
    ZZgt = overlap_poisson({V},{F},bc_fun,rhs_fun);
    av(s) = avgedge([VA;VB],[FA;FB+size(VA,1)]);
    avgt(s) = avgedge(V,F);
    errord(s) = max(abs(phi([VA;VB])-[ZZd{1};ZZd{2}]));
    errorn(s) = max(abs(phi([VA;VB])-[ZZn{1};ZZn{2}]));
    errorgt(s) = max(abs(phi(V)-ZZgt{1}));
    plot(log(av),log(errord),log(av),log(errorn),...
        log(avgt),log(errorgt),'LineWidth',3)
    legend('dirichlet','naive','groundtruth')
    title('Poisson equation convergence (nonsymmetric)')
    xlabel('log h')
    ylabel('log max error')
    axis equal
    drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - Silvia
