function [Aeq,II,outer_bb,inner_bb] = overlap_constraints(VV,FF,varargin)
  % OVERLAP_CONSTRAINTS
  %
  % [Aeq,II,outer_bb,inner_bb,C] = overlap_constraints(VV,FF,varargin)
  %
  % Inputs:
  %   VV  k-long list of #VV{i} by dim lists of vertex positions
  %   FF  k-long list of #FF{i} by dim+1 lists of element indices into VV{i}
  %   Optional:
  %     'Sparsification' followed by one of:
  %       'none'  generate constraints for every vertex in every mesh it lies in
  %       'first'  generate csonstraints for every vertex in the first mesh it
  %         lies in
  %       {'max-cover'}  **approximate** finding a maximum cover of vertices with
  %         the constraint that all vertices lying in at least one mesh are
  %         "covered" 
  %     'BoundaryFacets'  followed by precomputed k-long list of #OO{i} by dim
  %        list of boundary facet indices 
  %     'Method'  followed by one of:
  %       'naive','lsq'
  %       {'dirichlet'}
  %     'Dilate'  followed by number of dilation iterations on inner boundary
  %       {0}
  % Outputs:
  %   Aeq  #Aeq by sum(#VV{i}) equality constraints matrix
  %   II  k-long list of #VV{i} by k matrices so that II(i)(v,j) indicates
  %     whether vertex v of mesh i is inside mesh j.
  %   outer_bb   
  %   inner_bb   
  %   %C  k-long list of #VV{i} lists revealing the number of constraints
  %   %  "generated" by each vertex
  %

  k = numel(VV);
  assert(k == numel(FF));

  % default values
  method = 'dirichlet';
  sparsification = 'max-cover';
  dilate = 0;
  OO = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'BoundaryFacets','Dilate','Method','Sparsification'}, ...
    {'OO','dilate','method','sparsification'});
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
  dim = size(VV{1},2);

  IE = cell(k,1);
  for i = 1:k
    IE{i} = sparse(size(VV{i},1),k);
  end
  %V = cell2mat(VV{:;
  %ind_all = 1:size(V,1);
  for i = 1:k
    %ind_i = (n(i)+1):n(i+1);
    %ind_not_i = setddiff(ind_all,ind_i);
    %V_not_i = V(ind_not_i,:);
    others = [1:i-1 i+1:k]';
    V_not_i = cell2mat(arrayfun(@(j) VV{j},others,'UniformOutput',false));
    n_not_i = cumsum([0 cellfun(@(V) size(V,1),{VV{others}})]);
    switch dim
    case 1
      IEF = in_element(VV{i},FF{i},V_not_i);
      [IEI,IEJ] = find(IEF);
      IEi = zeros(size(V_not_i,1),1);
      IEi(IEI) = IEJ;
    otherwise
      IEi = in_element_aabb(VV{i},FF{i},V_not_i);
    end
    for j = 1:numel(others)
      IE{others(j)}(:,i) = IEi(n_not_i(j)+1:n_not_i(j+1));
    end
    % A vertex is considered to be inside its own mesh
    IE{i}(:,i) = 1;
  end


  % Collect running number of vertices in preceeding meshes
  n = cumsum([0 cellfun(@(V) size(V,1),{VV{:}})]);
  % Collect boundary edges of each mesh
  if isempty(OO)
    OO = cellfun(@(F) boundary_faces(F),FF,'UniformOutput',false);
  end
  bb = cellfun(@(O) unique(O),OO,'UniformOutput',false);
  %bbin = cell(k,k);
  %for i = 1:k
  %  others = [1:i-1 i+1:k]';
  %  % Boundary vertices of all other components
  %  Vb_not_i = cell2mat(arrayfun(@(j) VV{j}(bb{j},:),others,'UniformOutput',false));
  %  nb = cumsum([0 cellfun(@(b) numel(b),{bb{others}})]);
  %  if isempty(Vb_not_i)
  %    in_i = [];
  %  else
  %    switch dim
  %    case 1
  %      in_i = Vb_not_i > min(VV{i}(OO{i})) & Vb_not_i < max(VV{i}(OO{i}));
  %    otherwise
  %      % determine wether vertices of other mesh are inside this mesh
  %      %tic;
  %      in_i =  abs(winding_number(VV{i},OO{i},Vb_not_i))>0.5;
  %      %fprintf('winding_number (%d): %g secs\n',i,toc);
  %    end
  %  end
  %  % extract "insideness" for boundary vertices
  %  for j = 1:numel(others)
  %    bbin{others(j),i} = in_i(nb(j)+1:nb(j+1));
  %  end
  %  % a boundary vertex js not consjder to be inside its own mesh
  %  bbin{i,i} = false(numel(bb{i}),1);
  %end

  % truly outer boundary
  outer_bb = cell(k,1);
  % part of boundary inside any other shape
  inner_bb = cell(k,1);
  for i = 1:k
    bi = bb{i};
    inner = false(numel(bi),1);
    for j = [1:i-1 i+1:k]
      %i_in_j = abs(winding_number(VV{j},OO{j},VV{i}(bi,:)))>0.5;
      %inner = inner | bbin{i,j};
      inner = inner | IE{i}(bb{i},j);
    end
    others = [1:i-1 i+1:k]';
    outer_bb{i} = bi(~inner);
    inner_bb{i} = bi(inner);
    if dilate > 0
      Bi = sparse(inner_bb{i},1,1,size(VV{i},1),1);
      Ai = adjacency_matrix(FF{i}) + speye(size(VV{i},1));
      % include neighbors
      for iter = 1:dilate
        Bi = (Ai*Bi);
      end
      inner_bb{i} = find(Bi);
    end
  end

  II = cell(k,1);
  for i = 1:k
    II{i} = sparse(size(VV{i},1),k);
    II{i}(:,i) = 1;
  end

  Aeq = {};
  if k == 1
    Aeq = [];
    II = {true(size(VV{1},1),k)};
    C = {};
    return;
  end

  switch sparsification
  case 'first'
  % C{i}(v) reveals number of constraints that vertex v in mesh i has generated.
    C = cell(k,1);
    for i = 1:k
      C{i} = zeros(size(VV{i},1),1);
    end
  end

  for i = 1:k
    for j = [1:i-1 i+1:k]
      %[IC,C1,C2,C3] = in_element(VV{i},FF{i},VV{j}, ...
      %  'First',true,'Quiet',true,'Method','spatial-hash');

      % To-do: We should call this once for all non-i vertices so that the aabb
      % tree is not rebuilt k-1 times.
      %I = in_element_aabb(VV{i},FF{i},VV{j});
      %assert(all((I) == (IE{j}(:,i))))
      I = IE{j}(:,i);

      switch sparsification
      case 'first'
        I(C{j} > 0) = 0;
      end
      [I,~,J] = find(I);
      % We need to compute II (for all vertices, but we don't need B for all
      % vertices if doing 'dirichlet')
      II{j}(I,i) = 1;
      switch method
      case 'dirichlet'
        K = ismember(I,inner_bb{j});
        I = I(K);
        J = J(K);
      end
      switch size(VV{i},2)
      case 1
        B = barycentric_coordinates( ...
          VV{j}(I,:), ...
          VV{i}(FF{i}(J,1),:),VV{i}(FF{i}(J,2),:));
      case 2
        B = barycentric_coordinates( ...
          VV{j}(I,:), ...
          VV{i}(FF{i}(J,1),:),VV{i}(FF{i}(J,2),:),VV{i}(FF{i}(J,3),:));
      case 3
        B = barycentric_coordinates( ...
          VV{j}(I,:), ...
          VV{i}(FF{i}(J,1),:), ...
          VV{i}(FF{i}(J,2),:), ...
          VV{i}(FF{i}(J,3),:), ...
          VV{i}(FF{i}(J,4),:));
      end

      switch sparsification
      case 'first'
        C{j}(I) = C{j}(I) +1;
      end

      Aeqij =  ...
        sparse( ...
          repmat(1:numel(I),size(B,2),1)', ...
          FF{i}(J,:)+n(i), ...
          B, ...
          numel(I),n(end)) + ...
        sparse((1:numel(I))',I+n(j),-1,numel(I),n(end));
      if size(Aeqij,1) >0
        Aeq = {Aeq{:},Aeqij};
      end
    end
  end
  Aeq = cell2mat(Aeq');

  switch sparsification
  case 'max-cover'
    if true
      % This is not really a "max cover". A vertex may lie inside more than one
      % element from another mesh. Each generates a row in the constraint matrix
      % Aeq. We are going to keep only one such constraint per vertex. 
      Aeq_count = sum(Aeq>0,1);
      [I,J,V] = find(Aeq);
      score = max(sparse(I,J,Aeq_count(J),size(Aeq,1),size(Aeq,2)),[],2);
      [Im,Jm,Vm] = find(Aeq==-1);
      %[valid,Ir] = minnz(sparse(Im,Jm,score(Im),size(Aeq,1),size(Aeq,2)));
      inv_score = 1./score;
      [valid,Ir] = max(sparse(Im,Jm,inv_score(Im),size(Aeq,1),size(Aeq,2)),[],1);
      Aeq = Aeq(Ir(valid~=0),:);
    else
      Aeq_count = sum(Aeq~=0,1);
      remove = [];
      % For each column (vertices) of Aeq
      for i=1:size(Aeq,2)
        % Find all rows where this vertex is the "interpolated" one
        r = find(Aeq(:,i)==-1);
        if length(r)>1
          maxs = [];
          % Consider each constraint
          for j=1:length(r)
            % all other vertices involved in this constraint
            interp = find(Aeq(r(j),:));
            % don't consider the interpolated vertex itself, but actually do...
            maxs(j) = 0;
            % consider all vertices involved in constraints on interpolated vertex
            for s=1:length(interp)
              % Find all the rows that interp(s) appear in
              %srows = find(Aeq(:,interp(s))~=0);
              maxs(j) = max(maxs(j),Aeq_count(interp(s)));
            end
          end
          %r 
          %maxs
          %error
          % Sort in ascending order
          [~,ind] = sort(maxs);
          % Remove all constraints "on" this vertex except the one with the smallest
          % max-constrained other vertex.
          remove = [remove;r(ind(2:end))];
          
          %         s = find((Aeq(:,i)~=0).*(Aeq(:,i)~=-1));
          %         if length(s)>50
          %           remove = [remove; s(ceil(rand(1,floor(length(s)/10))*length(s)))];
          %         end
          %           Aeq(s(12:end),:)=[];
        end
      end
      Aeq(remove,:) = [];
    end
  end
  

  %switch sparsification
  %case 'first'
  %  remove = [];
  %  for i=1:size(Aeq,2)
  %    r = find(Aeq(:,i)==-1);
  %    if length(r)>1
  %      Aeq(r(2:end),:)=[];
  %      remove = [remove;r(2:end)];
  %    end
  %    s = find((Aeq(:,i)~=0).*(Aeq(:,i)~=-1));
  %    if length(s)>0
  %      disp('lol')
  %      k = [k,length(s)];
  %      if length(s)>40
  %        ind = [ind,i];
  %      end
  %      %                     Aeq(s(12:end),:)=[];
  %    end
  %  end
  %  Aeq(remove,:) = [];
  %end

end
