function [mV,mF,VV,FF] = qslim_components(V,F,t)
  % QSLIM_COMPONENTS
  %
  % [mV,mF,VV,FF] = qslim_components(V,F)
  %
  % Inputs:
  %   V  #V by 3 list of surface mesh vertex positions
  %   F  #F by 3 list of surface mesh face indices into V
  %   t  target number of faces {0 for all collapse}
  % Outputs:
  %   mV  #mV by 3 list of output surface mesh vertex positsions
  %   mF  #mF by 3 list of output surface mesh face indices into SV
  %

  C = connected_components(F); C = C(F(:,1));
  k = max(C);
  VV = cell(k,1);
  FF = cell(k,1);
  mV = cell2mat(VV);
  m = size(F,1);
  f = t/m;
  ms = zeros(k,1);
  for i = 1:k
    Fi = F(C==i,:);
    [Vi,I] = remove_unreferenced(V,Fi);
    Fi = I(Fi);
    VV{i} = Vi;
    FF{i} = Fi;
    ms(i) = size(Fi,1);
  end
  p = 0.5;
  %ms = m-ms;
  fs = (ms/sum(ms)).^p/sum((ms/sum(ms)).^p);
  for i = 1:k
    % Everybody does their part
    ti = ceil(fs(i)*t);
    [~,VV{i},FF{i}] = qslim(VV{i},FF{i},ti);
  end
  mV = cell2mat(VV);
  n = cumsum([0 cellfun(@(V) size(V,1),{VV{:}})]);
  mF = cell2mat(arrayfun(@(i) n(i)+FF{i},1:numel(VV),'UniformOutput',false)');
end

