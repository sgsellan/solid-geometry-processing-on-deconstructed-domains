function [V,T,F,SV,SF] = spherical_shell(s,R0,R1);
  % SPHERICAL_SHELL Generate the 3D equivalent of an annulus.
  %
  % [V,F] = spherical_shell(s,R0,R1)
  %
  % Inputs:
  %   s  number of vertices on inner shell
  %   R0  radius of inner shell
  %   R1  radius of outer shell
  % Outputs:
  %   V  #V by 3 list of vertix positions


  method = 'quartet';
  switch method
  case 'quartet'
    path_to_quartet = '/Users/ajx/Dropbox/quartet/quartet_release';
    quartet_output_file = 'quartet-output.tet';
    dx = min(1/sqrt(s),0.1*R1);
    obj_file = sprintf('quartet-sphere-%g-%g.obj',R0,R1);
    if ~exist(obj_file,'file')
      [SV,SF] = subdivided_sphere(7);
      SF = [fliplr(SF);size(SV,1)+SF];
      SV = [SV*R0;SV*R1];
      writeOBJ(obj_file,SV,SF);
    end
    cmd = sprintf('%s %s %0.17f %s',path_to_quartet,obj_file,dx,quartet_output_file);
    [status,result] = system(cmd);
    if status ~= 0
      fprintf('%s\n',cmd);
      error(result);
    end
    [V,T] = readTET(quartet_output_file);
    F = boundary_faces(T);
  case 'tetgen'
    [IV,IF] = lloyd_sphere(s);
    IV = IV*R0;
    IF = fliplr(IF);
    [OV,OF] = lloyd_sphere(ceil(s*(R1/R0)^2));
    OV = OV*R1;
    SV = [OV;IV];
    SF = [OF;size(OV,1)+IF];
    [V,T,F] = tetgen( ...
      SV,SF, ...
      'Flags',sprintf('-q1.2a%0.17f',avgedge(SV,SF)^3/(6*sqrt(2))),'Holes',[0 0 0]);
    % 24 comes from experimental  testing. It doesn't even seem to be enough...
    %[avgedge(SV,SF)^3/(6*sqrt(2)) median(volume(V,T))]
    %clf;
    %hold on;
    %tsurf(OF,OV,'FaceAlpha',0.5,'EdgeAlpha',0.5);
    %tsurf(IF,IV,'FaceAlpha',0.5,'EdgeAlpha',0.5);
    %hold off;
    %mean(doublearea(IV,IF))
    %mean(doublearea(OV,OF))
    %axis equal;
  case 'shells'
    h = 0.73*2*pi*R0/sqrt(s);
    rs = linspace(R0,R1,ceil((R1-R0)/h)+1);
    SV = cell(numel(rs),1);
    SF = cell(numel(rs),1);
    n = 0;
    max_iters = 10;
    [IV,IF] = lloyd_sphere(ceil(s),'MaxIters',max_iters);
    for ri = 1:numel(rs)
      r = rs(ri);
      ni = ceil(s * (r/R0)^2);
      [IV,IF] = lloyd_sphere(ni, ...
        'InitialGuess',[IV;randsphere(ni-size(IV,1))],'MaxIters',max_iters);
      SV{ri} = IV*r;
      if ri == 1 || ri == numel(rs)
        SF{ri} = n+IF;
      end
      n = n+size(IV,1);
    end
    SF = cell2mat(SF);
    SV = cell2mat(SV);
    [V,T,F] = tetgen( ...
      SV,SF, ...
      'Flags','-q','Holes',[0 0 0]);
  case 'cvt'
    n = ceil(sqrt(s)^3);
    U = rand(n,1).^(1/3);
    V = randsphere(n) .* (R0+U.*(R1-R0));
    T = delaunay(V);

    M = massmatrix(V,T);
    max_iters = 5;
    for iter = 1:max_iters
      A = adjacency_matrix(T);
      A = A*M;
      %A = bsxfun(@rdivide,A,sum(A,2));
      A = spdiags (1./sum (A,2), 0, size(A,1), size(A,1)) * A ;
      V_prev = V;
      V = A*V;
      % subtract off center of mass  (needed for small n)
      V = bsxfun(@minus,V,diag(M)'*V./sum(diag(M)));
      BC = barycenter(V,T);
      inner = normrow(BC)<=R0;
      V(T(inner),:) = normalizerow(V(T(inner),:))*R0;
      outer = normrow(BC)>=R1;
      V(T(outer),:) = normalizerow(V(T(outer),:))*R1;
      T = delaunay(V);
      M = massmatrix(V,T);
      er = trace((V-V_prev)'*M*(V-V_prev));
      %if vis
      %set(t,'Vertices',V,'Faces',F, ...
      %  ... 'CData',full(diag(M)));
      %  'CData',full(sum(adjacency_matrix(F),2)),'FaceLighting','phong','FaceColor','interp');
      %axis equal;
      %drawnow;
      %title(sprintf('%g',er));
      %end
      if er < 1e-07
        break;
      end
    end

    BC = barycenter(V,T);
    T = T(normrow(BC)>R0,:);
    [V,I] = remove_unreferenced(V,T);
    T = I(T);
    F = boundary_faces(T);
    SV = [];
    SF = [];
  end

end
