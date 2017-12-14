function PpS = plus(P,S, method)
%
%  PLUS: Add a Polyhedron and a vector or Polyhedron. 
%  ===================================================
%  
%  
%  SYNTAX
%  ------
%     
%      Q = P + S
%      Q = plus(P,S,method)
%      Q = P + x
%      Q = plus(P, x)
%    
%  
%  DESCRIPTION
%  -----------
%     Compute the Minkowski sum of P  and S, or P  and x. 
%                          P+S = { x+y  |  x in P, y in S}
%     or 
%                              P+y = { x+y  | x in P }
%  
%  
%  INPUT
%  -----
%     
%        
%          P      Polyhedron in any format                 
%                 Class: Polyhedron                        
%          S      Polyhedron in any format or a point      
%                 given as column real vector of the       
%                 length P.Dim                             
%                 Class: Polyhedron or double              
%          method String selecting the projection method   
%                 to be used. Can be 'vrep', 'fourier', or 
%                 'mplp'.                                  
%                 Class: string                            
%                   
%  
%  
%  OUTPUT
%  ------
%     
%        
%          Q Polyhedron Q = P+S.                      
%            Class: Polyhedron                        
%              
%  
%  
%  SEE ALSO
%  --------
%     minus
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
 
 
if nargin < 3
	% Default projection method. Empty=let Polyhedron/projection() decide.
	method = [];
end

if ~isa(P, 'Polyhedron')
	% v+P (issue #59)
	
	% note that "S" is for sure a Polyhedron object, or otherwise this
	% method would not even be called by matlab
	Sold = S;
	S = P;
	P = Sold;
end

% deal with arrays
if numel(P)>1
    if isa(S,'Polyhedron')
        PpS(numel(P)*numel(S),1) = Polyhedron;
        % match each polyhedron together
        for i=1:numel(P)
            for j=1:numel(S)
                PpS(numel(S)*(i-1)+j) = plus(P(i),S(j));
            end
        end
    else
        PpS(size(P)) = Polyhedron;
        for i=1:numel(P)
           PpS(i) = plus(P(i),S);
        end
    end
    return;
end

type = class(S);

switch type
    case 'Polyhedron'
        
        if numel(S)>1
            error('Only one polyhedron S is allowed.');
        end
                
        % if P is empty array
        if builtin('isempty',P)
            PpS = Polyhedron(S);
            return
        end        
        
        % Adding empty polyhedra
        if isEmptySet(P)
            PpS = Polyhedron(S);
            return
        end
        if isEmptySet(S)
            PpS = Polyhedron(P);
            return
        end
        
        if P.Dim ~= S.Dim,
            error('P and S must be the same dimension');
        end
        
        if P.hasVRep && S.hasVRep
            %% P is V-rep and S is V-rep
            Vp = P.V; Vs = S.V;
            Rp = P.R; Rs = S.R;
            
			if size(Vp, 1)==0
				V = Vs;
			elseif size(Vs, 1)==0
				V = Vp;
			else
				V = zeros(size(Vp,1)*size(Vs,1), P.Dim);
				for i = 1:size(Vs, 1)
					V((i-1)*size(Vp, 1)+1:i*size(Vp, 1), :) = bsxfun(@plus, Vp, Vs(i, :));
				end
			end
			
			if size(Rp, 1)==0
				R = Rs;
			elseif size(Rs, 1)==0
				R = Rp;
			else
				R = zeros(size(Rp,1)*size(Rs,1), P.Dim);
				for i = 1:size(Rs, 1)
					R((i-1)*size(Rp, 1)+1:i*size(Rp, 1), :) = bsxfun(@plus, Rp, Rs(i, :));
				end
			end
            
            PpS = Polyhedron('V',V,'R',R);
            
        elseif P.hasHRep && S.hasHRep
            %% P is H-rep and S is H-rep
            
            % Form sum
            Z  = zeros(size(P.A,1), size(S.A,2));
            Ze = zeros(size(P.Ae,1),size(S.Ae,2));
            tmp = Polyhedron('H', [S.A -S.A S.b;Z P.A P.b], 'He', [S.Ae -S.Ae S.be;Ze P.Ae P.be]);
            PpS = projection(tmp, 1:S.Dim, method);
            
        elseif P.hasHRep && S.hasVRep
            % P is H-rep and S is V-rep
            %%% TODO : Do this more efficiently
            
            % Do this the silly way now, by converting and adding
            P.minVRep();
            PpS = plus(P, S, method);
            
        elseif P.hasVRep && S.hasHRep
            % P is V-rep and S is H-rep
            
            %%% TODO : Do this more efficiently
            
            % Do this the silly way now, but converting and adding
            S.minVRep();
            PpS = plus(P, S, method);
        end
        
        % add function handle
        %PpS.Func = P.Func;

        
    otherwise
        
        validate_realvector(S);
        if size(S, 2)~=1 || numel(S)~=P.Dim
            error('The point must be a %dx1 vector', P.Dim);
        end
        x = S(:);
        
        % Get the data for P
        PV  = P.V_int;
        PR  = P.R_int;
        PH  = P.H_int;
        PHe = P.He_int;
        
        % V-rep - just shift the vertices (rays are not affected)
        if size(PV,1) > 0
            PV = PV + repmat(x', size(PV,1), 1);
        end
        
        % H-rep - A*z <= b, Ae*z == be
        %  ==> Ay <= b + A*x,  Ae*y == be + Ae*x
        if size(PH,1) > 0
            PH(:,end) = PH(:,end) + PH(:,1:end-1)*x;
        end
        if size(PHe,1) > 0
            PHe(:,end) = PHe(:,end) + PHe(:,1:end-1)*x;
        end
        
        % Minkowski sum with a vector does not affect irredundancy
		if isempty(PV) && isempty(PR) && isempty(PHe) && ~isempty(P.H_int)
			% use a specific syntax for faster execution
			PpS = Polyhedron(PH(:, 1:end-1), PH(:, end));
			PpS.irredundantHRep = P.irredundantHRep;
			PpS.irredundantVRep = P.irredundantVRep;
		else
			PpS = Polyhedron('V',PV,'R',PR,'H',PH,'He',PHe,...
				'irredundantVRep',P.irredundantVRep,'irredundantHRep',P.irredundantHRep);
		end
        %PpS.Func = P.Func;

end


end
