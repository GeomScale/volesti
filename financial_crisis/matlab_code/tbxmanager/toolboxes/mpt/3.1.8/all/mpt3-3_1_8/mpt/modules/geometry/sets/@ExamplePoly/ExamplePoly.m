classdef ExamplePoly
  % Create various random polyhedra
  
  methods(Static=true)
    %%
    function P = randHrep(varargin)
      ip = inputParser;
      ip.KeepUnmatched = false;
      ip.addParamValue('d',   2, @isnumeric);
      ip.addParamValue('n',   10, @isnumeric);
      ip.addParamValue('nr',  0, @isnumeric);
      ip.addParamValue('ne',  0, @isnumeric);
      ip.parse(varargin{:});
      p = ip.Results;

      % Random Hrep containing the origin the n inequalities, d-ne dimensions
      P = Polyhedron('H', [randn(p.n,p.d-p.ne) ones(p.n,1)]);
      
      % If an unbounded polyhedron is requested, then add a cone
      if p.nr > 0
        R = Polyhedron('R', randn(p.nr,p.d-p.ne));
        P = P + R;
        % Convert to inequality form
        P.minHRep();
      end
      n = size(P.A,1);
      
      % Lift and rotate into d dimensions
      A  = [P.A randn(n,p.ne)];
      Ae = [zeros(p.ne,p.d-p.ne) eye(p.ne)];
      [U,S,V] = svd(randn(p.d));
      
      P = Polyhedron('H', [A*U P.b], 'He', [Ae*U zeros(p.ne,1)]);
    end
    
    %%
    function P = randVrep(varargin)
      ip = inputParser;
      ip.KeepUnmatched = false;
      ip.addParamValue('d',   2, @isnumeric);
      ip.addParamValue('n',   10, @isnumeric);
      ip.addParamValue('nr',  0, @isnumeric);
      ip.addParamValue('ne',  0, @isnumeric);
      ip.parse(varargin{:});
      p = ip.Results;
      
      % Random Vrep containing the origin in Rd with n vertices and nr rays
      P = Polyhedron('V', randn(p.n,p.d-p.ne), 'R', randn(p.nr,p.d-p.ne));
      
      % If ne is non-zero, then rotate P up into d-dimensional space
      [U,S,V] = svd(randn(p.d));
      P = P.affineMap(U(:,1:p.d-p.ne));
    end
    
    %%
    function P = randZono(varargin)
      ip = inputParser;
      ip.KeepUnmatched = false;
      ip.addParamValue('d',   2, @isnumeric);
      ip.addParamValue('n',   5, @isnumeric);
      ip.parse(varargin{:});
      p = ip.Results;
      
      P = Polyhedron;
      for i = 1:p.n
        v = randn(1,p.d);
        P = P + Polyhedron('V',[v;-v]);
        P.minVRep();
      end
    end
    
    %% 
    function P = poly3d_sin(varargin)
      ip = inputParser;
      ip.KeepUnmatched = false;
      ip.addParamValue('d',   2, @isnumeric);
      ip.addParamValue('n',   13, @isnumeric);
      ip.parse(varargin{:});
      p = ip.Results;
      n = p.n; d = p.d;
      
      % H-rep polytope whose normals are sinusoids
      th = linspace(0,2*pi,n+1)'; th(end) = [];
      N = [sin(th) cos(th)];
      A = zeros(0,d);
      for i=1:d
        for j=i+1:d
          Anew = zeros(n,d);
          Anew(:,i) = N(:,1); Anew(:,j) = N(:,2);
          A = [A;Anew];
        end
      end
      
      P = Polyhedron('H',[A ones(size(A,1),1)]);
    end
    
  end
  
end
