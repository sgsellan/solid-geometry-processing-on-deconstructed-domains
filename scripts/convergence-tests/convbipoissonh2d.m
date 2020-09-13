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
% 'strict':
%
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
    ZZd = overlap_bipoissonh(VV,FF,bc_fun,rhs_fun,norm_der,'Method','dirichlet');
    ZZn = overlap_bipoissonh(VV,FF,bc_fun,rhs_fun,norm_der,'Method','naive');
    ZZs = overlap_bipoissonh(VV,FF,bc_fun,rhs_fun,norm_der,'Method','strict');
    ZZgt = overlap_bipoissonh({V},{F},bc_fun,rhs_fun,norm_der);
    av(s) = avgedge([VA;VB],[FA;FB+size(VA,1)]);
    avgt(s) = avgedge(V,F);
    errord(s) = max(abs(phi([VA;VB])-[ZZd{1};ZZd{2}]));
    errorn(s) = max(abs(phi([VA;VB])-[ZZn{1};ZZn{2}]));
    errors(s) = max(abs(phi([VA;VB])-[ZZs{1};ZZs{2}]));
    errorgt(s) = max(abs(phi(V)-ZZgt{1}));
    plot(log(av),log(errord),log(av),log(errorn),log(av),log(errors),...
        log(avgt),log(errorgt),'LineWidth',3)
    legend('dirichlet+laplacian','naive+laplacian','dirichlet+fullhessian','groundtruth')
    title('Bilaplace equation convergence')
    xlabel('log h')
    ylabel('log max error')
    axis equal
    drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - Silvia