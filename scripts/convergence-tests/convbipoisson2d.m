%%%%%%%%%%%%%%-2D BIPOISSON EQUATION CONVERGENCE TEST-%%%%%%%%%%%%%%%%%%%%%
%   This script is designed to exhaustively test convergence of the 2D
%bipoisson equation with Dirichlet and Neumann boundary conditions
%                         lap lap u = f,
%               u|_{\partial\Omega} = g,
%           du/dn|_{\partial\Omega} = h,
%where f,g,h are known functions (rhs_fun, bc_fun and norm_der respectively)
%and an element discretization of \Omega is not known. Instead, we use element
%discretizations of intersecting domains \Omega_i where
%                    \cup\Omega_i = \Omega.
%   We will test two different 'intersection' constraints: 'naive' and
%'dirichlet', which are in essence equality constraints on the intersection
%or the boundary of the intersection (inner boundaries) of the domains,
%respectively. Mathematically, this amounts to for all j,i=1,...,k,
% 'naive':    
%            u_i|_{\Omega_i\cap\Omega_j}=u_j|_{\Omega_i\cap\Omega_j}
%       \lap u_i|_{\Omega_i\cap\Omega_j}=\lap u_j|_{\Omega_i\cap\Omega_j}
% 'dirichlet':    
%  u_i|_{\partial(\Omega_i\cap\Omega_j)}=u_j|_{\partial(\Omega_i\cap\Omega_j)}
%\lap u_i|_{\partial(\Omega_i\cap\Omega_j)}=\lap u_j|_{\partial(\Omega_i\cap\Omega_j)}
%
%   In most example we also compute the numerical 'ground-truth', where we
%solve on a discretization of the whole \Omega, for reference.
%
%
%
%%%%%%%%%%%%%%%%%%%% BILAPLACE EQUATION ON AN ANNULUS %%%%%%%%%%%%%%%%%%%%%
%   We consider \Omega to be an annulus of radii 1 and 2, which we can see
%as the intersecting union of \Omega_1 and \Omega_2, annuli of radii 1 and
%1.8 and 1.2 and 2, respectively. In this case, we solve the equation
%                             bilap u = 0
% subject to
%             u|_{R=1} = 1                 u|_{R=2} = 0
%         du/dn|_{R=1} = 0             du/dn|_{R=2} = 0
%%%%%% Finding analytical solution with separation of variables %%%%%%%%%%%
clear
hold off
figure(1)
a = 1;
b = 2;
Z = @(r) (-4*a^2*b^2*log(b)^2 + 2.*b^2*log(b)*(a^2 - r.^2 + 2*a^2*log(r)) + ...
  2.*a^2*log(a)*(-b^2 + r.^2 + 2.*b^2*log(b) - 2.*b^2*log(r)) - (a^2 - ...
  b^2)*(b^2 - r.^2 + 2*r.^2.*log(r)))/((a^2 - b^2)^2 - ...
  4.*a^2*b^2*log(a)^2 + 8.*a^2*b^2*log(a)*log(b) - 4.*a^2*b^2*log(b)^2);
phi = @(V) Z(normrow(V));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Loop to find and plot errors %%%%%%%%%%%%%%%%%%%%%%%%%%%
bc_fun = @(V) normrow(V)<1.5;
rhs_fun = @(V) zeros(size(V,1),1);
norm_der = rhs_fun;
for s=1:6
    [VA,FA] = annulus(2^(s+3),2,'R',1.2);
    [VB,FB] = annulus(2^(s+3),1.8,'R',1);
    [V,F] = annulus(2^(s+3),2,'R',1);
    VV = {VA,VB};
    FF = {FA,FB};
    ZZd = overlap_bipoisson(VV,FF,bc_fun,rhs_fun,norm_der,'Method','dirichlet');
    ZZn = overlap_bipoisson(VV,FF,bc_fun,rhs_fun,norm_der,'Method','naive');
    ZZgt = overlap_bipoisson({V},{F},bc_fun,rhs_fun,norm_der);
    av(s) = avgedge([VA;VB],[FA;FB+size(VA,1)]);
    avgt(s) = avgedge(V,F);
    errord(s) = max(abs(phi([VA;VB])-[ZZd{1};ZZd{2}]));
    errorn(s) = max(abs(phi([VA;VB])-[ZZn{1};ZZn{2}]));
    errorgt(s) = max(abs(phi(V)-ZZgt{1}));
    plot(log(av),log(errord),log(av),log(errorn),...
        log(av),log(errorgt),'LineWidth',3)
    legend('dirichlet','naive','groundtruth')
    title('Bilaplace equation convergence')
    xlabel('log h')
    ylabel('log max error')
    axis equal
    drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% BIPOISSON EQUATION ON AN ANNULUS %%%%%%%%%%%%%%%%%%%%%
%   We consider \Omega to be an annulus of radii 1 and 2, which we can see
%as the intersecting union of \Omega_1 and \Omega_2, annuli of radii 1 and
%1.8 and 1.2 and 2, respectively. In this case, we will sample many
%different functions, their bilaplacians and their normal derivatives on
%the annulus.
clear
hold off
phi = @(V) (V(:,1).^4).*(V(:,2).^4);
bc_fun = phi;
rhs_fun = @(V) 24.*((V(:,1).^4)+(V(:,2).^4)+12.*(V(:,1).^2).*(V(:,2).^2));
theta = @(V) cart2pol(V(:,1),V(:,2));
norm_der = @(V) 2.*(-4.*(sqrt(V(:,1).^2+V(:,2).^2)<1.5).*(cos(theta(V)).^4).*(sin(theta(V)).^4)+...
    512.*(sqrt(V(:,1).^2+V(:,2).^2)>1.5).*(cos(theta(V)).^4).*(sin(theta(V)).^4));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi2 = @(V) sin(V(:,1));
bc_fun2 = phi2;
rhs_fun2 = phi2;
norm_der2 = @(V) (-(sqrt(V(:,1).^2+V(:,2).^2)<1.5).*(cos(cos(theta(V))).*cos(theta(V)))+...
    (sqrt(V(:,1).^2+V(:,2).^2)>1.5).*(cos(2.*cos(theta(V))).*cos(theta(V))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi3 = @(V) sin(V(:,1)).*sin(V(:,2));
bc_fun3 = phi3;
rhs_fun3 = @(V) 4.*sin(V(:,1)).*sin(V(:,2));
norm_der3 = @(V) ((-1).*(sqrt(V(:,1).^2+V(:,2).^2)<1.5).*((cos(cos(theta(V))).*sin(sin(theta(V))).*cos(theta(V)))+...
    sin(cos(theta(V))).*cos(sin(theta(V))).*sin(theta(V)))+(sqrt(V(:,1).^2+V(:,2).^2)>1.5).*((cos(2.*cos(theta(V))).*sin(2.*sin(theta(V))).*cos(theta(V)))+...
    sin(2.*cos(theta(V))).*cos(2.*sin(theta(V))).*sin(theta(V))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Loop to find and plot errors %%%%%%%%%%%%%%%%%%%%%%%%%%%
for s=1:6
    [VA,FA] = annulus(2^(s+3),2,'R',1.2);
    [VB,FB] = annulus(2^(s+3),1.8,'R',1);
    [V,F] = annulus(2^(s+3),2,'R',1);
    VV = {VA,VB};
    FF = {FA,FB};
    ZZd = overlap_bipoisson(VV,FF,bc_fun,rhs_fun,norm_der,'Method','dirichlet');
    ZZn = overlap_bipoisson(VV,FF,bc_fun,rhs_fun,norm_der,'Method','naive');
    ZZgt = overlap_bipoisson({V},{F},bc_fun,rhs_fun,norm_der);
    av(s) = avgedge([VA;VB],[FA;FB+size(VA,1)]);
    avgt(s) = avgedge(V,F);
    errord(s) = max(abs(phi([VA;VB])-[ZZd{1};ZZd{2}]));
    errorn(s) = max(abs(phi([VA;VB])-[ZZn{1};ZZn{2}]));
    errorgt(s) = max(abs(phi(V)-ZZgt{1}));
    ZZd2 = overlap_bipoisson(VV,FF,bc_fun2,rhs_fun2,norm_der2,'Method','dirichlet');
    ZZn2 = overlap_bipoisson(VV,FF,bc_fun2,rhs_fun2,norm_der2,'Method','naive');
    ZZgt2 = overlap_bipoisson({V},{F},bc_fun2,rhs_fun2,norm_der2);
    errord2(s) = max(abs(phi2([VA;VB])-[ZZd2{1};ZZd2{2}]));
    errorn2(s) = max(abs(phi2([VA;VB])-[ZZn2{1};ZZn2{2}]));
    errorgt2(s) = max(abs(phi2(V)-ZZgt2{1}));
    ZZd3 = overlap_bipoisson(VV,FF,bc_fun3,rhs_fun3,norm_der3,'Method','dirichlet');
    ZZn3 = overlap_bipoisson(VV,FF,bc_fun3,rhs_fun3,norm_der3,'Method','naive');
    ZZgt3 = overlap_bipoisson({V},{F},bc_fun3,rhs_fun3,norm_der3);
    errord3(s) = max(abs(phi3([VA;VB])-[ZZd3{1};ZZd3{2}]));
    errorn3(s) = max(abs(phi3([VA;VB])-[ZZn3{1};ZZn3{2}]));
    errorgt3(s) = max(abs(phi3(V)-ZZgt3{1}));
end
figure(2)
plot(log(av),log(errord),log(av),log(errorn),...
        log(av),log(errorgt),'LineWidth',3)
legend('dirichlet','naive','groundtruth')
title('Bipoisson equation convergence test 1')
xlabel('log h')
ylabel('log max error')
axis equal
figure(3)
plot(log(av),log(errord2),log(av),log(errorn2),...
        log(av),log(errorgt2),'LineWidth',3)
legend('dirichlet','naive','groundtruth')
title('Bipoisson equation convergence test 2')
xlabel('log h')
ylabel('log max error')
axis equal
figure(4)
plot(log(av),log(errord3),log(av),log(errorn3),...
        log(av),log(errorgt3),'LineWidth',3)
legend('dirichlet','naive','groundtruth')
title('Bipoisson equation convergence test 3')
xlabel('log h')
ylabel('log max error')
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - Silvia