function [ZZ,data,UU] = overlap_poisson(VV,FF,bc_fun,rhs_fun,varargin)
% OVERLAP_POISSON  Solve a Poisson equation on a domain composed by taking the
% union of overlapping meshes (VV,FF). (Dirichlet) boundary conditions and
% right-hand side values are provided via a pointwise evaluation funtion.
%
%_
%
% ZZ = overlap_poisson(VV,FF,bc_fun,rhs_fun)
% ZZ = ...
%   overlap_poisson(VV,FF,bc_fun,rhs_fun,'ParameterName',ParameterValue, ...)
%
% Inputs:
%   VV  k-long list of #VV{i} by dim lists of vertex positions
%   FF  k-long list of #FF{i} by ss lists of element indinces into VV{i}
%   bc_fun  boundary conditions funtion so that bc_fun(x) returns the
%     boundary condition value at position x. Returning nan indicates that
%     this position should not receive an explicit value constraint (i.e., it
%     should receive natural boundary conditions instead)
%   rhs_fun  Poisson equation right-hand side funtion so that rhs_fun(x)
%     returns the pointwise right-hand side value at position x (not yet
%     integrated).
%   Optional:
%     'Alpha'  followed by penalty weight (only used for 'Method','penalty')
%       {1000}
%     'ConstrainAuxiliary'  followed by wether or not to constraint auxilliary
%       variables durign mixed FEM solve {true}
%     'Data'  (see output)
%     'Dilate'  see overlap_constraints.m
%     'Domain'  csg function (not used)
%     'FixedIndices'  followed by #fixed_indices list of indices into
%       cat(1,VV{:}) of fixed vertices (does not work for 'Method','corefine')
%       {[]}
%     'FixedValues'  followed by #fixed-indices list of values (does not work
%       for 'Method','corefine') {[]}
%     'K'  Order of k-Poisson equation to solve
%     'Lambda' solve M+λL instead. **Warning** the mass matrix used here is in
%       real units. You might see better numerics by scaling down input models
%       so that mass matrix is closer to [0,1] range. In which case, you'll
%       have to also scale down your Lambda value correspondingly {inf}
%      'LowerBounds'  followed by  #V
%     'LowerBound'/'UpperBound'  followed by lower and upper bounds respectively
%     'Method' followed by one of:
%       'corefine'  combine the overlapping shapes into a single mesh
%       'dirichlet' constrain only boundary vertices
%       'penalty'  constrain all vertcies (like 'naive') but with a soft
%         penalty objective
%       'naive'  constrain all vertices of all meshes to linear combination of
%         corners of triangles they fall in of other meshes
%     'Solver'  followed by one of:
%       'ipqpcommon'  matlab's qp solver. This seems pretty good for large
%          inputs
%       'min_quad_with_fixed'  this may trigger QR factorization and then become
%         unbearably slow
%       'min_quad_with_fixed-force_Aeq_li'  this is probably _not_ ok for
%         'Method','naive'. But probably OK for 'Method','dirichlet'
%       'quadprog'  use either matlab or mosek qp solver
%       By default if K==1 then {'ipqpcommon'} and otherwise
%         {'min_quad_with_fixed'}
%     'Sparsification'  see overlap_constraints.m
%     
% Outputs:
%   ZZ  k-long list of #VV{i} solutions at vertices of each mesh
%   data  struct containing precomputed data
%    .GG  k-long list of #FF{i}*dim by #VV{i} gradient matrices
%    .AA  k-long list of #FF{i} by #FF{i} diagonal _adjusted_ area matrices
%    .MM  k-long list of #VV{i} by #VV{i} barycenter _adjusted_ mass matrices
%    .LL k-long list of #VV{i} by #VV{i} barycenter _adjusted_ -Laplace matrices
%    .corefine  corefine data
%    .Aeq  #Aeq by #VV overlap constraints matrix
%    .aux_Aeq  #aux_Aeq by #VV auxilliary constraints (when K>1)
%

  function [b,bc] = filter_bc(V,b)
    bc = bc_fun(V(b,:));
    keep = ~isnan(bc);
    bc = bc(keep,:);
    b = b(keep);
  end

  k = numel(VV);
  dim = size(VV{1},2);
  ss = size(FF{1},2);
  assert(k == numel(FF));

  % default values
  fixed_indices = [];
  fixed_values = [];
  lb = [];
  ub = [];
  t = 1; %test
  u=[];
  interp = 0;
  interpverts = [];
  h = 0;
  msbk = 0;
  boundedweights = 0;
  waveeq = 'false';
  method = 'dirichlet';
  csg = 'union';
  constrain_auxiliary = true;
  sparsification = 'max-cover';
  penalty_alpha = 1000;
  lambda = inf;
  data = [];
  dilate = 0;
  solver = [];
  vertices = [];
  K = 1;
  samples = 0;
  mc = 0;
  domain_fun = @csgunion;
  bchack = 0;
  collect_timings = true;
  % union
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
      { ...
        'Alpha','ConstrainAuxiliary','Data','Dilate','Domain', ...
        'FixedIndices','FixedValues','K',...
        'Lambda','LowerBound','Method','Solver','Sparsification', ...
        'UpperBound', ...
      ...  %undocumented, please ignore:
      'Samples','MC','Domain_fun','CSG',...
      'BoundedWeights','Vertices','WaveEq','H','U','Interp','InterpData','MSBK','T',...
      'BChack'}, ...
      {'penalty_alpha','constrain_auxiliary','data','dilate','domain_fun', ...
       'fixed_indices','fixed_values','K',...
       'lambda','lb','method','solver','sparsification','ub', ...
      ...  %undocumented, please ignore:
      'samples','mc','domain_fun','csg',...
      'boundedweights','vertices','waveeq','h','u','interp','interpdata','msbk','t',...
      'bchack'});
  v = 1;
  while v <= numel(varargin)
      param_name = varargin{v};
      if isKey(params_to_variables,param_name)
          assert(v+1<=numel(varargin));
          v = v+1;
          % Trick: use feval on anonymous funtion to use assignin to this workspace
          feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
      else
          error('Unsupported parameter: %s',varargin{v});
      end
      v=v+1;
  end

  if collect_timings
    fprintf('#T: %d\n',size(cat(1,FF{:}),1));
    fprintf('K: %d\n',numel(FF));
  end

  assert(numel(fixed_indices) == size(fixed_values,1));

  if isempty(solver)
    switch method
    %case 'dirichlet'
    %  % This does seem faster for big meshes, but slower than ipqpcommon for
    %  % smaller meshes. But I'm still not 100% sure that data.Aeq is truly, always
    %  % full rank...
    %  if size(data.Q,1) > 50000
    %    solver = 'min_quad_with_fixed-force_Aeq_li';
    %  else
    %   solver = 'ipqpcommon';
    %  end
    otherwise
      if isempty(lb) && isempty(ub) 
        if K == 1
          solver = 'ipqpcommon';
        else
          solver = 'min_quad_with_fixed';
        end
      else
        solver = 'quadprog';
      end
    end
  end

  if isempty(data)
    has_data = false;
    data = struct('Method',method);
  else
    has_data = true;
    assert(isfield(data,'Method'));
  end

  % Collect running number of vertices in preceeding meshes
  n = cumsum([0 cellfun(@(V) size(V,1),{VV{:}})]);
  if ~has_data
    % Collect boundary edges of each mesh
    O = cellfun(@(F) boundary_faces(F),FF,'UniformOutput',false);
  end

  switch method
  case 'corefine'
    assert(K == 1);
    assert(~has_data);
    % Gather all vertices in one mesh
    OV = cell2mat(VV(:));
    E = [];
    b = [];
    % boundary edges indexing the joined mesh
    for i = 1:k
        E = [E;n(i)+O{i}];
        b = [b;n(i)+unique(O{i})];
    end
    switch dim
    case 2
      % Create a single mesh containing all vertices (and possibly vertices at
      % boundary edge intersections)
      [V,F] = triangle(OV,E,[],'Flags','');
      % Only keep the triangles inside the original shapes
    case 3
      [SV,SE] = selfintersect(OV,E,'StitchAll',true);
      %% mesh_boolean will throw away unreferenced points...
      % [SV,SE] = mesh_boolean(OV,E,[],[],'union');
      % clf;
      % hold on;
      % tsurf(SE,SV,'Facealpha',0.8,'FaceColor',blue);
      % scatter3(SV(:,1),SV(:,2),SV(:,3),'.');
      % hold off;
      % axis equal;
      % axis([-1 1 -1 0 -1 1])
      % drawnow;
      [V,F] = tetgen(SV,SE,'Flags','');
      %medit(V,F);
    end
    W = winding_number(OV,E,barycenter(V,F));
    F = F(abs(W)>0.5,:);
    % Get the actual boundary of this shape
    b = unique(boundary_faces(F));
    % Build the Laplace and mass matrix
    Q = -cotmatrix(V,F);
    M = massmatrix(V,F);
    l = [2*M*rhs_fun(V)];
    [b,bc] = filter_bc(V,b);
    data.mqwf = struct();
    data.mqwf.force_Aeq_li = true;
    Z = min_quad_with_fixed(Q,l,b,bc,[],[],data.mqwf);
    % Distribute solution: if there are duplicates in OV then they won't all
    % show up in V so need to find mapping.
    I = knnsearch(V,OV,'K',1);
    ZZ = cell(k,1);
    for i = 1:k
      ZZ{i} = Z(I(n(i)+1:n(i+1)),:);
    end
    data.corefine.V = V;
    data.corefine.F = F;
    return;
  end

  % Precompute data used for solve
  if ~has_data
    tic;
    [data.Aeq,II,bb,inner_bb] = ...
      overlap_constraints(VV,FF, ...
        'Dilate',dilate, ...
        'BoundaryFacets',O, ...
        'Method',method,'Sparsification',sparsification);
    time_constraints = toc;
    if collect_timings
      fprintf('constraints: %g secs\n',time_constraints);
    end
    tic;
    [data.LL,data.MM,data.AA,data.GG] = ...
      overlap_operators(VV,FF,II,'MC',mc,'Samples',samples,'Domain',domain_fun);
    time_operators = toc;
    if collect_timings
      fprintf('operators: %g secs\n',time_operators);
    end

    data.Q = blkdiag(data.LL{:});
    data.L = data.Q; %for K=2
    Beq = zeros(size(data.Aeq,1),1);

    linear = [];
    data.b = [];
    % loop over each subdomain and gather rhs/linear term
    for i = 1:k
      data.b = [data.b;bb{i}+n(i)];
      linear = [linear;2*data.MM{i}*rhs_fun(VV{i})];
    end

    switch K
    case 2
      [data.Q,linear,data.Aeq,Beq] = overlap_bilaplace_system( ...
        data.LL,data.MM,linear,data.Aeq,Beq, ...
        'ConstrainAuxiliary',constrain_auxiliary,'InnerBB',inner_bb);
    end
    % These became unused ...
    data.aux_Aeq = [];
    aux_Beq = [];

    if ~isinf(lambda)
      assert(K == 1);
      data.Q = lambda*data.Q + blkdiag(data.MM{:});
    end

    switch method
    case 'penalty'
      assert(K == 1);
      % convert constraints into energy term
      data.Q = data.Q +   data.Aeq'*data.Aeq*penalty_alpha;
      linear = linear - 2*data.Aeq'*Beq*penalty_alpha;
      data.Aeq = [];
      Beq = [];
    end
    % Append constraints on auxilliary dofs
    data.Aeq = [data.Aeq;data.aux_Aeq];
    Beq = [Beq;aux_Beq];
  end

  V = cat(1,VV{:});
  [data.b,data.bc] = filter_bc(V,data.b);
  % enforce fixed values constraints
  data.b = [data.b(:);fixed_indices(:)];
  data.bc = [data.bc;fixed_values];

  % dummy bounds for now
  LB = -inf(size(data.Q,1),1);
  UB = inf(size(data.Q,1),1);
  if ~isempty(lb)
    LB(1:size(V,1)) = lb;
  end
  if ~isempty(ub)
    UB(1:size(V,1)) = ub;
  end

  % SOLVE The resulting quadratic programming problem:
  %
  % min           Z'*data.Q*Z + Z'*linear
  %  Z
  % subject to    data.Aeq Z = Beq 
  %               Z(data.b) = data.bc
  %               lb ≤ Z ≤ ub
  %
  % This should not contain any logic except selecting and calling the solver of
  % choice.
  tic;
  switch solver
  case {'min_quad_with_fixed','min_quad_with_fixed-force_Aeq_li'}
    if ~isfield(data,'mqwf')
      data.mqwf = struct();
    end
    switch solver
    case 'min_quad_with_fixed-force_Aeq_li'
      assert( ...
        ~isfield(data,'mqwf') || ...
        ~isfield(data.mqwf,'force_Aeq_li') || ...
        data.mqwf.force_Aeq_li);
      data.mqwf.force_Aeq_li = true;
    end
    Z = min_quad_with_fixed(data.Q,linear,...
      data.b,data.bc,data.Aeq,Beq,data.mqwf);    
  case 'ipqpcommon'
    Z = ipqpcommon_wrapper( ...
        data.Q,0.5*linear,[],[], ...
        [data.Aeq; ...
        sparse(1:numel(data.b),data.b,1,numel(data.b),size(data.Q,1))], ...
        [Beq;data.bc], ...
        LB,UB, zeros(0,1));
  case 'quadprog'
    param = default_quadprog_param;
      %param.Diagnostics = 'on';
      %param.Display = 'iter';
    Z = quadprog( ...
      data.Q,0.5*linear,[],[], ...
      [data.Aeq; ...
        sparse(1:numel(data.b),data.b,1,numel(data.b),size(data.Q,1))], ...
        [Beq;data.bc], ...
        LB,UB,[],param);
  end
  %ZU = reshape(Z,[],2)
  %ZU([3 4 5 11 12],:)
  time_solver = toc;
  if collect_timings
    fprintf('solve: %g secs\n',time_solver);
  end

  % Place solution into cell array
  switch K
  case 2
    % only keep primary variables
    %assert(n(end) == size(Z,1)/2);
    %U = Z(n(end)+(1:n(end)));
    Z = Z(1:n(end));
  end
  ZZ = cell(k,1);
  %UU = cell(k,1);
  for i = 1:k
    ZZ{i} = Z(n(i)+1:n(i+1),:);
    %if K>1
    %  UU{i} = U(n(i)+1:n(i+1),:);
    %end
  end

end
