function [V,T,F,C,n,Z] = combine_meshes(VV,TT,FF,ZZ,ks)
  % COMBINE MESHES combinatorially
  %
  % [V,T,F,C,n,Z] = combine_meshes(VV,TT,FF,ZZ,ks)
  % 
  % Inputs:
  %   VV  k-long list of #VV{i} by dim lists of vertex positions
  %   TT  k-long list of #TT{i} by 4 lists of face indices into VV{i}
  %   FF  k-long list of #FF{i} by 3 lists of face indices into VV{i}
  %   ZZ  k-long list of #VV{i} by dim list of vertex values
  % Outputs
  %   V  n(end) by dim list of vertex positions
  %   T  #T by 3 list of face indices into V
  %   F  #F by 3 list of face indices into V
  %   C  #C list of indices into 1:numel(ks)
  %   n  k+1-long list of cumulative sum of #VV{i} values
  %   Z  #V list of vertex values
  %
  % See also: overlap_constraints, overlap_poisson
  %
  if nargin <= 4 
    ks = 1:numel(VV);
  end
  V = cell2mat(VV(ks));
  n = cumsum([0 cellfun(@(V) size(V,1),{VV{ks}})]);
  F = cell2mat(arrayfun(@(i) n(i)+FF{ks(i)},1:numel(ks),'UniformOutput',false)');
  T = cell2mat(arrayfun(@(i) n(i)+TT{ks(i)},1:numel(ks),'UniformOutput',false)');
  C = cell2mat(arrayfun(@(i) repmat(i,size(VV{ks(i)},1),1),1:numel(ks),'UniformOutput',false)');
  if nargin >=4
    Z = cell2mat(ZZ(ks));
  end
end
