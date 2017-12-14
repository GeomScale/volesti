function Pnew = affineMap(H,T,method)
%
%  AFFINEMAP: Compute the affine map of the Polyhedron. 
%  =====================================================
%  
%  
%  SYNTAX
%  ------
%     
%      Q = P.affineMap(T)
%      Q = P.affineMap(T,method)
%      Q = affineMap(P,T,method)
%    
%  
%  DESCRIPTION
%  -----------
%     Computes an affine map P |->> Q  of polyhedron P to polyhedron Q based on the
%  transformation matrix T. The polyhedron Q is given by 
%                               n                    d          
%                     Q = { yinR   |  y = Tx, xin P R  }     (1)
%                                                               
%     The matrix T must be real with n rows and d columns. 
%    
%     - If n<d then this operation is referred to as projection. 
%     - If n=d then this operation is referred to as rotation/skew. 
%     - If n>d then this operation is referred to as lifting. 
%  
%  
%  INPUT
%  -----
%     
%        
%          P      Polyhedron in any format.                
%                 Class: Polyhedron                        
%          T      Transformation matrix.                   
%                 Class: double                            
%          method Specific method to use in projection     
%                 operation. Allowed methods are "vrep",   
%                 "fourier", and "mplp". For details type  
%                 "help Polyhedron/projection".            
%                 Class: string                            
%                   
%  
%  
%  OUTPUT
%  ------
%     
%        
%          Q Polyhedron representing the affine map   
%            in H- or V-representation.               
%            Class: Polyhedron                        
%              
%  
%  
%  SEE ALSO
%  --------
%     projection,  affineHull
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
%   mailto:colin.jones@epfl.ch 
%  
%  

%  LICENSE
%  -------
%    
%    This program is free software; you can redistribute it and/or modify it under
%  the terms of the GNU General Public License as published by the Free Software
%  Foundation; either version 2.1 of the License, or (at your option) any later
%  version.
%    This program is distributed in the hope that it will be useful, but WITHOUT
%  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
%  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%    You should have received a copy of the GNU General Public License along with
%  this library; if not, write to the  Free Software Foundation, Inc.,  59 Temple
%  Place, Suite 330,  Boston, MA 02111-1307 USA
%  -------------------------------------------------------------------------------
%    
%      This document was translated from LaTeX by HeVeA (1).
%  ---------------------------------------
%    
%    
%   (1) http://hevea.inria.fr/index.html
 
 
global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

if nargin < 3,
    method = [];
end

% deal with arrays
if numel(H)>1
    Pnew(size(H)) = Polyhedron;
    for i=1:numel(H)
        Pnew(i) = H(i).affineMap(T,method);
    end
    return;
end

% check dimension of T
validate_realmatrix(T);
if size(T,2)~=H.Dim
    error('The matrix representing the affine map must have %d number of columns.',H.Dim);
end

if H.isEmptySet,
	Pnew = Polyhedron.emptySet(size(T, 1));
    return;

elseif norm(T) <= MPTOPTIONS.abs_tol
	% Special case: zero map. Return a singleton (the origin) of
	% appropriate dimension: T*P = { z | z=0 } where the dimension of "z"
	% is equal to the number of rows of "T"
	new_dim = size(T, 1);
	V = zeros(new_dim, 1);
	He = [eye(new_dim), zeros(new_dim, 1)];
	Pnew = Polyhedron('V', V', 'He', He);
	return
end

if H.hasVRep
  Pnew = Polyhedron('V',H.V*T', 'R', H.R*T');
else
  
  if size(H.Ae,1) > 0
      
      if size(H.A, 1)==0
          % edge case: no inequalities
          if size(T, 1)==size(T, 2) && abs(det(T)) > MPTOPTIONS.abs_tol
              % simple solution if "T" is invertible
              Pnew = Polyhedron('Ae', H.Ae*inv(T), 'be', H.be);
          else
              % lower-dimensional mapping, solution as suggester by Magnus
              % Nilsson
              nT = size(T,1);
              p = H.Ae\H.be; % particular solution to Ae*d=be
              Ae = zeros(nT); % a bit ugly way to avoid empty matrices for Ae when it should be zero matrices.
              Ae = [null([T*null(H.Ae)]')';Ae]; % Concatenate with the zero matrix...
              Ae = Ae(1:nT, :);                 % ...and remove superfluous rows.
              be = Ae*T*p;
              Pnew = Polyhedron('Ae', Ae, 'be', be);
          end
          
      else
          % Remove equalities
          opt = Opt('P',H,'pF',T','f',zeros(size(T,2),1));
          opt = opt.eliminateEquations;
          
          % Compute the mapping for the new representation of the polyhedron
          Pnew = Polyhedron('H', [opt.A opt.b]).affineMap(opt.pF',method) + T*opt.recover.th(:,end);
      end
      return
  end

  % Compute permutation of T s.t. P*y = [T1 T2; T3 T4]*[xr;xn] with rank T1 = rank T and Q*x=[xr;xn]
  [L,U,p,q] = lu(sparse(T),'vector');
  r = rank(T,MPTOPTIONS.abs_tol);
  pr = p(1:r); pn = p(r+1:end);
  qr = q(1:r); qn = q(r+1:end);
 
  rk = rank(full(U(:,1:r)));
  if rk~=r
      % if invertibility is not achieved try reduced echelon elimination
      [~,jb]=rref(T,MPTOPTIONS.abs_tol);
      [L,U,p] = lu(sparse(T(:,jb)),'vector');
      % update column selection
      q = [jb, setdiff(1:H.Dim,jb)];
      pr = p(1:r); pn = p(r+1:end);
      qr = q(1:r); qn = q(r+1:end);
  end
  
  beta = T(pr,qr) \ eye(r);
  
  % Build polyhedron who's projection is the full-dimensional portion of the mapping
  A = H.A;
  Ptmp = Polyhedron(A(:,q)*[beta -beta*T(pr,qn); zeros(length(qn),r) eye(length(qn))], H.b);
  ptmp = Ptmp.projection(1:r,method);
  
  % Don't know what format will come back from the projection
  if ptmp.hasVRep
    V = zeros(size(ptmp.V,1),size(T,1));
    R = zeros(size(ptmp.R,1),size(T,1));
    V(:,pr) = ptmp.V;
    V(:,pn) = ptmp.V*(T(pn,qr)*beta)';
    R(:,pr) = ptmp.R;
    R(:,pn) = ptmp.R*(T(pn,qr)*beta)';
    
    Pnew = Polyhedron('V', V, 'R', R);
  else
    A  = zeros(size(ptmp.A,1),size(T,1));
    Ae = zeros(length(pn),size(T,1));
    A(:,pr) = ptmp.A;
    b = ptmp.b;
    Ae(:,pr) = T(pn,qr)*beta;
    Ae(:,pn) = -eye(length(pn));
    be = zeros(size(Ae,1),1);
    
    Pnew = Polyhedron('H', [A b], 'He', [Ae be]);
  end
  
%   % See the note "Affine Mapping of Polyhedra with Empty Interiors" by C.N.
%   % Jones for details on the algorithm.
%   
%   E = H.Ae;
%   c = H.be;
%   G = H.A;
%   g = H.b;
%   
%   % Compute permutation matrices via LU factorization and form partition
%   %  P[T;E]Q = [A B; C D], where rank A == rank [A B; C D]
%   [L,U,P,Q,R] = lu(sparse([T;E]));
%   r = rank([T;E], MPTOPTIONS.abs_tol);
%   tmp = P*[T;E]*Q;
%   [m,n] = size([T;E]);
%   mz = size(T,1); me = size(E,1);
%   A = tmp(1:r,1:r);     B = tmp(1:r,r+1:end);
%   C = tmp(r+1:end,1:r); D = tmp(r+1:end,r+1:end);
%   
%   if rank(A) ~= r, error('Could not find full-rank basis'); end
%   iA = inv(A);
%   
%   % Compute affine hull of the projection
%   F = [-C*iA eye(m-r)] * P * [eye(mz);zeros(me,mz)];
%   f = [-C*iA eye(m-r)] * P * [zeros(mz,1);c];
%   
%   % Take full-rank subset of rows of F
%   nF = []; nf = [];
%   for i=1:size(F,1)
%     if rank([nF nf;F(i,:) f(i)]) > rank([nF nf])
%       nF = [nF;F(i,:)];
%       nf = [nf;f(i)];
%     end
%   end
%   F = nF;
%   f = nf;
%   
%   % Form matrices of full-dimensional polyhedron to project
%   Gbar = G*Q*[iA*[eye(r) zeros(r,m-r)]*P*[eye(mz);zeros(me,mz)];zeros(n-r,mz)];
%   Lbar = G*Q*[iA*B;eye(n-r)];
%   gbar = g - G*Q*[iA*[eye(r) zeros(r,m-r)]*P*[zeros(mz,1);c];zeros(n-r,1)];
%   
%   % Compute partition [Fh Ft] = F*Pbar s.t. Fh is full-rank
%   if ~isempty(F)
%     [L,U,Qbar,Pbar,R] = lu(sparse(F));
%     rF = rank(F, MPTOPTIONS.abs_tol);
%     Fh = F*Pbar(:,1:rF);
%     Ft = F*Pbar(:,rF+1:end);
%     
%     if rank(Fh) ~= rank(F), error('Could not find full-rank basis for F'); end
%     
%     % Compute projection of full-dim polyhedron to full-dim polyhedron
%     Pfull = Polyhedron([Gbar*Pbar*[-inv(Fh)*Ft;eye(size(Ft,2))] Lbar], ...
%       gbar - Gbar*Pbar*[inv(Fh)*f;zeros(size(Pbar,2)-rF,1)]);
%     
%     pfull = Pfull.projection(1:size(Ft,2), method);
%     pfull.minHRep();
%     L = pfull.A; l = pfull.b;
%     
%     % Compute final result
%     Pnew = Polyhedron('He',[F f],'H',[L*[zeros(size(L,2),size(Pbar,1)-size(L,2)) eye(size(L,2))]*inv(Pbar) l]);
%   else
%     
%     % Compute final result
%     Pnew = projection(Polyhedron('H',[Gbar Lbar gbar]), 1:size(Gbar,2), method);
%   end
end
