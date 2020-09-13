function [LL,MM,AA,GG] = overlap_operators(VV,FF,II,varargin)
  % OVERLAP_OPERATORS Build discrete differential geometry operators for
  % overlapping domains.
  % 
  % [LL,MM,AA,GG] = overlap_operators(VV,FF,II,varargin)
  %
  % Inputs:
  %   VV  k-long list of #VV{i} by dim lists of vertex positions
  %   FF  k-long list of #FF{i} by ss lists of element indinces into VV{i}
  %   II  k-long list of #VV{i} by k matrices so that II(i)(v,j) indicates
  %     whether vertex v of mesh i is inside mesh j.
  % Outputs:
  %   LL k-long list of #VV{i} by #VV{i} barycenter _adjusted_ -Laplace matrices
  %   MM  k-long list of #VV{i} by #VV{i} barycenter _adjusted_ mass matrices
  %   AA  k-long list of #FF{i} by #FF{i} diagonal _adjusted_ area matrices
  %   GG  k-long list of #FF{i}*dim by #VV{i} gradient matrices
  %
  % See also: overlap_constraints, overlap_poisson
  mc = 0;
  samples = 0;
  domain_fun = [];
  params_to_variables = containers.Map( ...
      { 'Samples','MC','Domain'}, ...
      {'samples','mc','domain_fun'});
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

  k = numel(VV);
  dim = size(VV{1},2);
  ss = size(FF{1},2);
  LL = cell(k,1);
  MM = cell(k,1);
  AA = cell(k,1);
  GG = cell(k,1);
  % Loop over each subdomain and gather adjusted area terms and then build
  % Laplacian for each subdomain
  for i = 1:k
    if samples==0
      Vcount = sum(II{i},2);
      Fcount = mean(Vcount(FF{i}),2);
    else
      if mc==0
        V = VV{i};
        assert(size(V,2)==3);
        F = FF{i};
        [QQ,W] = quadraturepoints(V,F,samples);
        Fcount = sparse(size(F,1),samples);
        for r=1:samples
          Fcountr = sparse(size(F,1),k);
          for j=setdiff(1:k,i)
            Fcountr(:,j) = abs(winding_number(VV{j},O{j},QQ{i}))>.5;
          end
          Fcount(:,r) = W(r)*(sum(Fcountr,2)+1);
        end
        Fcount = sum(Fcount,2);
      else
        alpha = samples;
        barysamples1 = rand(alpha,1).^3;
        V = VV{i};
        F = FF{i};
        barysamples2 = rand(alpha,1).*(1-barysamples1);
        barysamples4 = rand(alpha,1).*(1-barysamples1-barysamples2);
        barysamples3 = 1-barysamples2-barysamples1-barysamples4;
        barys = [barysamples1, barysamples2, barysamples3,barysamples4];
        Fcount = sparse(size(F,1),alpha);
        for r=1:alpha
          switch dim
          case 3
            Ftest = V(F(:,1),:).*barysamples1(r)+V(F(:,2),:).*barysamples2(r)+...
            V(F(:,3),:).*barysamples3(r)+V(F(:,4),:).*barysamples4(r);
          case 2
            Ftest = V(F(:,1),:).*barysamples1(r)+V(F(:,2),:).*barysamples2(r)+...
            V(F(:,3),:).*barysamples3(r);
          end
          Fcountr = sparse(size(F,1),k);
          for j=setdiff(1:k,i)
            Fcountr(:,j) = abs(winding_number(VV{j},O{j},Ftest))>.5;
          end
          Fcount(:,r) = sum(Fcountr,2)+1;
        end
        Fcount = mean(Fcount,2);
      end
    end
    GG{i} = grad(VV{i},FF{i});
    switch dim
    case 1
      A = abs(VV{i}(FF{i}(:,2))-VV{i}(FF{i}(:,1)));
    case 2
      A = doublearea(VV{i},FF{i})*0.5;
    case 3
      A = volume(VV{i},FF{i});
    end
    % Doesn't make a big difference
    %Fcount = max(Vcount(FF{i}),[],2);
    A = A./Fcount;
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % A = A.*ficare{i};
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %Vdomain = domain_fun(II{i});
    %Fdomain = mean(Vdomain(FF{i}),2);
    %A(Fdomain<0.5) = eps*A(Fdomain<0.5);
    AA{i} = A;
    % Barycentric mass matrix
    MM{i} = sparse(FF{i},FF{i},repmat(A,1,ss)/ss,size(VV{i},1),size(VV{i},1));
    LL{i} = GG{i}'*repdiag(diag(sparse(A)),dim)*GG{i};
  end

end

