%%%%%%%%%%%%%%-2D POISSON EQUATION CONVERGENCE TEST-%%%%%%%%%%%%%%%%%%%%%%%
%   This script is designed to exhaustively test convergence of the 2D poisson
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
%%%%%%%%%%%%%%%%%%%% LAPLACE EQUATION ON AN ANNULUS %%%%%%%%%%%%%%%%%%%%%%%
%   We consider \Omega to be an annulus of radii 1 and 2, which we can see
%as the intersecting union of \Omega_1 and \Omega_2, annuli of radii 1 and
%1.8 and 1.2 and 2, respectively. In this case, we solve the equation
%                             lap u = 0
% subject to
%             u|_{R=1} = 1                 u|_{R=2} = 0
%
%%%%%% Finding analytical solution with separation of variables %%%%%%%%%%%
clear
hold off
figure(1)
A=[1,log(1);1,log(2)];
b=[1;0];
C=inv(A)*b;
C1=C(1);
C2=C(2);
phi = @(V) C1+C2*log(normrow(V));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Loop to find and plot errors %%%%%%%%%%%%%%%%%%%%%%%%%%%
bc_fun = @(V) normrow(V)<1.5;
rhs_fun = @(V) zeros(size(V,1),1);
for s=1:7
    [VA,FA] = annulus(2^(s+4),2,'R',1.2);
    [VB,FB] = annulus(2^(s+4),1.8,'R',1);
    [V,F] = annulus(2^(s+4),2,'R',1);
    VV = {VA,VB};
    FF = {FA,FB};
    ZZd = overlap_poisson(VV,FF,bc_fun,rhs_fun,'Method','dirichlet');
    ZZn = overlap_poisson(VV,FF,bc_fun,rhs_fun,'Method','naive');
    %ZZgt = overlap_poisson({V},{F},bc_fun,rhs_fun);
    av(s) = avgedge([VA;VB],[FA;FB+size(VA,1)]);
    %avgt(s) = avgedge(V,F);
    errord(s) = max(abs(phi([VA;VB])-[ZZd{1};ZZd{2}]));
    errorn(s) = max(abs(phi([VA;VB])-[ZZn{1};ZZn{2}]));
    %errorgt(s) = max(abs(phi(V)-ZZgt{1}));
    plot(log(av),log(errord),log(av),log(errorn),...
        'LineWidth',3)
    legend('dirichlet','naive')
    title('Laplace equation convergence')
    xlabel('log h')
    ylabel('log max error')
    axis equal
    drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% POISSON EQUATION ON AN ANNULUS %%%%%%%%%%%%%%%%%%%%%%%
%   We consider \Omega to be an annulus of radii 1 and 2, which we can see
%as the intersecting union of \Omega_1 and \Omega_2, annuli of radii 1 and
%1.8 and 1.2 and 2, respectively. In this case, we solve the equation
%                             lap u = 1
%   subject to
%                       u|_{R=1} = u|_{R=2} = 0
%
%%%%%% Finding analytical solution with separation of variables %%%%%%%%%%%
clear
hold off
figure(2)
R0 = 1;
R1 = 2;
phir = @(r) (-1/4).*(-r.^2+(R0.^2-R1.^2)./(log(R0)-log(R1)).*log(r)+...
    (log(R0).*R1.^2-R0.^2.*log(R1))./(log(R0)-log(R1)));
phi = @(V) phir(normrow(V));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Loop to find and plot errors %%%%%%%%%%%%%%%%%%%%%%%%%%%
bc_fun = @(V) zeros(size(V,1),1);
rhs_fun = @(V) ones(size(V,1),1); 
for s=1:5
    [VA,FA] = annulus(2^(s+4),2,'R',1.2);
    [VB,FB] = annulus(2^(s+4),1.8,'R',1);
    [V,F] = annulus(2^(s+4),2,'R',1);
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
        log(av),log(errorgt),'LineWidth',3)
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
clear
hold off
figure(3)
phi = @(V) (V(:,1).^2).*(V(:,2).^2);
bc_fun = phi;
rhs_fun = @(V) (2.*(V(:,2).^2))+(2.*(V(:,1).^2)); %!!!!!!!!!!!!!!!!!!!!
for s=1:5
    [VA,FA] = annulus(2^(s+4),2,'R',1.2);
    [VB,FB] = annulus(2^(s+4),1.8,'R',1);
    [V,F] = annulus(2^(s+4),2,'R',1);
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
        log(av),log(errorgt),'LineWidth',3)
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