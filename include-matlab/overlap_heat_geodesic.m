function ZZ = overlap_heat_geodesic(VV,TT,source_fun,varargin)
  % OVERLAP_HEAT_GEODESIC
  %
  % ZZ = ...
  %   overlap_heat_geodesic(VV,TT,source_fun,lambda,'ParamName',ParamValue, ...)
  %
  % Inputs:
  %   VV  k-long list of #VV{i} by dim lists of vertex positions
  %   FF  k-long list of #FF{i} by ss lists of element indinces into VV{i}
  %   source_fun  function source_fun(V) returning whether position in each row
  %     should be consider a source.
  %   Optional:
  %     'Lambda'  followed by time step parameter of "Geodesics in Heat" (solve
  %       using (M+λL)
  %     'BoundaryConditions'  followed by one of the following strings
  %       {'average'}: 
  %       'dirichlet'  domain boundaries set to 0 when computing heat
  %         diffusion
  %       'neumann'  heat diffusion solved with implicit neumann conditions
  % Outputs:
  %   ZZ  k-long list of #VV{i} solutions at vertices of each mesh
  %
  
  % default values
  overlap_method = 'dirichlet';
  lambda = 0.1;
  bc_type = 'average';
  sparsification = 'max-cover';
  solver = 'min_quad_with_fixed-force_Aeq_li';
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'BoundaryConditions','Lambda','OverlapMethod','Solver','Sparsification'}, ...
    {'bc_type','lambda','overlap_method','solver','sparsification'});
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

  switch bc_type
  case {'dirichlet','average'}
    bc_fun = @(V) source_fun(V);
    rhs_fun = @(V) -source_fun(V);
    [UUd,opdata] = overlap_poisson( ...
    VV,TT,bc_fun,rhs_fun, ...
    'Method',overlap_method,'Lambda',lambda,'Solver',solver,'Sparsification',sparsification);
    Ud = cell2mat(UUd);
  end
  
  switch bc_type
  case {'neumann','average'}
    bc_fun = @(V) nan(size(V,1),1);
    rhs_fun = @(V) -source_fun(V);
    [UUn,opdata] = overlap_poisson( ...
      VV,TT,bc_fun,rhs_fun,'Method',overlap_method,'Solver',solver,'Lambda',lambda, ...
      'Sparsification',sparsification);
    Un = cell2mat(UUn);
  end
  
  switch bc_type
  case 'neumann'
    UU = UUn;
  case 'dirichlet'
    UU = UUd;
  case 'average'
    UU = arrayfun(@(i) (UUd{i}+UUn{i})/2,1:numel(UUd),'UniformOutput',false)';
  end
  U = cell2mat(UU);
  
  
  % remove this and use opdata
  k = numel(VV);
  DD = cell(k,1);
  opdata.GG = cell(k,1);
  for i = 1:k
    opdata.GG{i} = grad(VV{i},TT{i});
    DD{i} = div(VV{i},TT{i});
    DD{i} = -0.5*opdata.GG{i}'*repdiag(diag(sparse(opdata.AA{i})),size(VV{i},2));
  end
  G = blkdiag(opdata.GG{:});
  D = blkdiag(DD{:});
  X = cell2mat( ...
    arrayfun(@(i) ...
      reshape(normalizerow(reshape(opdata.GG{i}*UU{i},[],3)),[],1), ...
      1:k,'UniformOutput',false)');
  % Some G*U rows may have had zero length
  X(any(isinf(X) | isnan(X),2),:) = 0;
  div_X = D*X;
  
  tic;
  Z = min_quad_with_fixed( ...
    blkdiag(opdata.LL{:}) + 1e-5*blkdiag(opdata.MM{:}),-2*div_X, ...
    [],[],opdata.Aeq,zeros(size(opdata.Aeq,1),1), ...
    struct('force_Aeq_li',true));
    fprintf('  min_quad_with_fixed: %g secs\n',toc);
  ZZ = cell(k,1);
  n = cumsum([0 cellfun(@(V) size(V,1),{VV{:}})]);
  for i = 1:k
    ZZ{i} = Z(n(i)+1:n(i+1),:);
  end
end

