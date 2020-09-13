function [Q,linear,Aeq,Beq] = overlap_bilaplace_system(LL,MM,linear,Aeq,Beq,varargin)
  % [Q,linear,Aeq,Beq] = overlap_bilaplace_system(LL,MM,linear,Aeq,Beq,varargin)
  %
  % Inputs:
  %   LL  k-long list of #VV{i} by #VV{i} barycenter _adjusted_ -Laplace matrices
  %   MM  k-long list of #VV{i} by #VV{i} barycenter _adjusted_ mass matrices
  %   linear  sum(#VV{i}) by #linear list of current linear coefficients (rhs_fun)
  %   Aeq  #Aeq by sum(#VV{i}) current equality constraints matrix
  %   Beq  #Beq by #linear current equality constraints rhss
  %   Optional:
  %     'ConstrainAuxiliary'  followed by wether or not to constraint auxilliary
  %       variables durign mixed FEM solve {true}
  %     'InnerBB'  #k list of #inner_bb{i} list of inner boundary indices
  %       (needed for 'ConstrainAuxiliary',false)  
  % Outputs (assuming #:
  %   Q  #Q by #Q square matrix representing bilaplace system/energy
  %   linear  #Q by #linear linear coefficients 
  %   Aeq  #Aeq by #Q  constraints matrix (perhaps modified from input)
  %   Beq  #Beq by #linear  cosntraints rhs
  %

  v = 1;
  constrain_auxiliary = true;
  inner_bb = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'ConstrainAuxiliary','InnerBB'}, ...
    {'constrain_auxiliary','inner_bb'});

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

  k = numel(LL);

  %M = blkdiag(MM{:});
  M = sparse(0,0);
  for i = 1:k
    Mi = MM{i};
    M = blkdiag(M,Mi);
    if ~constrain_auxiliary
      % get natural boundary conditions instead of zero-neumann (see [Stein
      % et al. 2017])
      LL{i}(:,inner_bb{i}) = 0;
    end
  end

  % This is the legacy version. I'm leaving it here so that we can compute
  % the comparison figure where we don't constrain the auxiliarry variables.
  legacy = ~constrain_auxiliary;
  if legacy
    %% We'll solve for Z and U=∆Z via [Jacobson et al. 2010]
    %Q = blkdiag(sparse(size(M,1),size(M,1)),M);
    %linear = [linear;zeros(size(M,1),size(linear,2))];
    %aux_Aeq = [blkdiag(LL{:})' -M];
    %aux_Beq = zeros(size(aux_Aeq,1),size(linear,2));

    % This works but the result is not a convex energy so can't use quadprog
    L = blkdiag(LL{:});
    % The transpose is important when ~constrain_auxiliary
    Q = [sparse(size(M,1),size(M,1)) L;L' -M];
    linear = [linear;zeros(size(M,1),size(linear,2))];
    aux_Aeq = [];
    aux_Beq = [];

    % Constraints due to overlapping domains
    if constrain_auxiliary
      Aeq = repdiag(Aeq,2);
      Beq = [Beq;Beq];
    else
      % Only constrain primary dofs
      Aeq = [Aeq sparse(size(Aeq,1),size(Aeq,2))];
    end
  else
    L = blkdiag(LL{:});

    % we would like to solve:
    %
    % saddle [z u]'[0 L;L -M][z u]  subject to A z = 0, A u = 0
    %  z,u
    %
    % but we need it as an energy so we can call quadprog. So we introduce a
    % Lagrange multiplier and then factor out u:
    %
    % min ‖L z - A' γ‖²_(M⁻¹) subject to A z = 0 
    % z,γ
    %

    %Q = [L Aeq']'*(M\[L Aeq']);
    %% This works but has a nasty little parameter
    %%Ik = speye(size(Aeq,1),size(Aeq,1));
    %%Q = [ ...
    %%  L*(M\L)                 (L/M)*Aeq'; ...
    %%  Aeq*(M\L)  Ik*1e-5+Aeq*(M\(Aeq'))];
    %Aeq = [Aeq sparse(size(Aeq,1),size(Aeq,1))];
    %linear = [linear;zeros(size(Aeq,1),size(linear,2))];

    % Sigh. Why is life so hard? The M⁻¹ above can get really nasty for tet
    % meshes and then Q becomes poorly conditioned. So we transform this
    % problem with another auxiliary variable y:
    %
    %  min  ‖y‖² subject to A z = 0, L z + A' γ = √M y
    % z,γ,y
    % 
    nq = size(L,1);
    na = size(Aeq,1);
    Q = blkdiag(sparse(nq+na,nq+na),speye(nq,nq));
    linear = [linear;zeros(na+nq,size(linear,2))];
    Msqrt = diag(sqrt(diag(M)));
    Aeq = [ ...
      Aeq sparse(na,na+nq);
      [L Aeq'] -Msqrt];
    Beq = [Beq;zeros(nq,size(linear,2))];
  end


end
