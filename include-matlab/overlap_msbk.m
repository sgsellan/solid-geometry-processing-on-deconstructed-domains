function ZZ = overlap_msbk(VV,FF,t,b)
  % Multi-Scale Biharmonic Kernel (actually the "prebiharmonic kernel" at the
  % end of Section 5 of "Multiscale Biharmonic Kernels" [Rustamov 2011])
  %
  % ZZ = overlap_msbk(VV,FF,t,b)
  %
  % Inputs:
  %   VV  k-long list of #VV{i} by dim lists of vertex positions
  %   FF  k-long list of #FF{i} by ss lists of element indinces into VV{i}
  %   t  time parameter
  %   b  (I'm being lazy, this is a list of boundary vertex ids. This should be
  %     a bc_fun)
  % Outputs:
  %   ZZ  k-long list of #VV{i} solutions at vertices of each mesh

  collect_timings = true;

  OO = cellfun(@(F) boundary_faces(F),FF,'UniformOutput',false);
  [Aeq,II,outer_bb,inner_bb] = overlap_constraints(VV,FF,'BoundaryFacets',OO);

  if collect_timings
    fprintf('overlap_constraints: %g secs\n',timeit(@() overlap_constraints(VV,FF,'BoundaryFacets',OO)));
  end
  [LL,MM] = overlap_operators(VV,FF,II);
  M = blkdiag(MM{:});
  ne = size(Aeq,1);
  n = size(Aeq,2);
  [Q,linear,Aeq,Beq] =  ...
    overlap_bilaplace_system(LL,MM,zeros(n,1),Aeq,zeros(ne,1));
  lb = [zeros(n,1);-inf(n+ne,1)];
  ub = [ones(n,1);  inf(n+ne,1)];
  lb(b) = 1;
  ub(b) = 1;
  tic;
  tic;
  fr = quadprog( ...
    Q,linear, ...
    [diag(M)' zeros(1,n+ne)],t, ...
    Aeq,Beq, ...
    lb,ub, ...
    [], ...
    default_quadprog_param);
  if collect_timings
    fprintf('quadprog: %g\n',toc);
  end
  Z = fr(1:n);

  k = numel(VV);
  ZZ = cell(k,1);
  n = cumsum([0 cellfun(@(V) size(V,1),{VV{:}})]);
  for i = 1:k
    ZZ{i} = Z(n(i)+1:n(i+1),:);
  end
end
