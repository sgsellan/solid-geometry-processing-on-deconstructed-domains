function [VV,TT,FF,I] = tetgen_components(V,F,varargin)
  % TETGEN_COMPONENTS  Given a multi-component and overlapping (but otherwise
  % clean) surface triangle mesh, construct tet-meshes of each connected
  % component independently.
  %
  % [VV,TT,FF,II] = tetgen_components(V,F)
  %
  % Inputs:
  %   V  #V by 3 list of surface mesh vertex positions
  %   F  #F by 3 list of surface mesh face indices into V
  % Outputs:
  %   VV  #components list of tet-mesh vertex positions
  %   TT  #components list of tet-mesh tet indices into VV{i}
  %   FF  #components list of tet-mesh boundary face indices into FF{i}
  %   I  #V list of indices into Vcat := cell2mat(VV) so that Vcat(I,:) = V
  %
  % Example:
  %   [VV,TT,FF,I] = tetgen_components(OV,OF,'Flags','-q100');
  %   V = cell2mat(VV);
  %   n = cumsum([0 cellfun(@(V) size(V,1),{VV{:}})]);
  %   F = cell2mat(arrayfun(@(i) n(i)+FF{i},1:numel(VV),'UniformOutput',false)');
  %   subplot(2,1,1);
  %   % Render tetmesh boundary faces (using tetmesh vertices)
  %   tsurf(F,V);
  %   % Render original faces BUT using tetmesh vertices!
  %   subplot(2,1,2);
  %   tsurf(I(OF),V);
  %   

  % default values
  max_vol_factor = 1;
  flags = [];
  P = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Flags','MaxVolFactor','Points'}, ...
    {'flags','max_vol_factor','P'});
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

  C = connected_components(F); C = C(F(:,1));
  if isempty(flags)
    max_vol = max_vol_factor*avgedge(V,F)^3/(sqrt(2));
    flags = sprintf('-q1.2a%0.17f',max_vol);
  end
  k = max(C);
  I = 1:size(V,1);
  VV = cell(k,1);
  TT = cell(k,1);
  FF = cell(k,1);
  II = cell(k,1);
  JJ = cell(k,1);
  n = 0;
  for i = 1:k
    Fi = F(C==i,:);
    [Vi,II{i},JJ{i}] = remove_unreferenced(V,Fi);
    I(JJ{i}) = n+(1:size(Vi,1));
    Fi = II{i}(Fi);
    try
      [VV{i},TT{i},FF{i}] = tetgen([Vi;P],Fi, ...
        'Flags',flags);
       n = n+size(VV{i},1);
    catch err
      fprintf('Failed on %d\n',i);
      rethrow(err);
    end

  end
end
