function N = nan_if(L)
  % NAN_IF Create a matrix with nans at non-zero entries of input matrix
  %
  % N = nan_if(L)
  %
  % Inputs:
  %   L  m by n matrix of logicals
  % Outputs:
  %   N  m by n matrix so that N(i,j) is nan if L(i,j) is true (non-zero)
  %
  N = zeros(size(L));
  N(L) = nan;
end
