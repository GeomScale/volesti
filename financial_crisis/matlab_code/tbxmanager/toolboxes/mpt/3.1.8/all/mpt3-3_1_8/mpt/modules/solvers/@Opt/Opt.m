classdef Opt < handle & matlab.mixin.Copyable
%
%  OPT: Interface for solving optimization problems 
%  =================================================
%  
%  
%  SYNTAX
%  ------
%     
%      problem = Opt('H',H,'f',f,'pF',F,'A',A,'b',b,'Ae',Ae,...)
%      problem = Opt('M',M,'q',q,'Q',Q,'Ath',Ath,'bth',bth)
%      problem = Opt(Matrices)
%      problem = Opt(constraints, objective, theta, x)
%    
%  
%  DESCRIPTION
%  -----------
%     Encapsulates data and solutions for LP/QP/MPLP/MPQP/LCP/PLCP problems. The
%  data can be provided directly, or using YALMIP format or MPT2 format. The
%  following general MPQP format that is considered by Opt class: 
%                     1  T             T         T                               
%               min   - x Hx+(Ftheta+f) x + theta Ytheta + Ctheta + c        (1) 
%                     2                                                          
%              s.t.   Ax <= b + Btheta                                       (2) 
%                                                                                
%                     A x = b  + Etheta                                      (3) 
%                      e     e                                                   
%                                                                                
%                     A     theta = b                                        (4) 
%                      theta         theta                                       
%     which contains n  optimization variables x, d  number of parameters y,
%  minequality constrains, m_e  equality constraints and constraints on the
%  parameter theta. Opt class accepts also PLCP formulation of the form 
%                                w - Mz = q + Qtheta        (5) 
%                                             w >= 0        (6) 
%                                             z >= 0        (7) 
%                                             T                 
%                                            w z = 0        (8) 
%                                                               
%     where w  and z  are the optimization variables and theta  is the parameter.
%  The output format of the optimization problem depends on the input data. For
%  instance, non-parametric LCP can be created by providing M, q  data and the
%  parametric LPC by supplying M, q, and Q. At the time of construction, the
%  problem data are validated and checked for proper dimensions. For all the
%  methods involved in Opt class see methods('Opt').
%  
%  INPUT
%  -----
%     
%        
%          H           Quadratic part of the objective          
%                      function.                                
%                      Class: double                            
%                      Default: []                              
%          f           Linear part of the objective function.   
%                      Class: double                            
%          pF          Linear part of the objective function    
%                      for parameters.                          
%                      Class: double                            
%                      Default: []                              
%          Y           Quadratic part of the objective function 
%                      for parameters.                          
%                      Class: double                            
%                      Default: 0                               
%          C           Linear part of the objective function    
%                      for parameters.                          
%                      Class: double                            
%                      Default: 0                               
%          c           Constant term in the objective function  
%                      Class: double                            
%                      Default: 0                               
%          A           Linear part of the inequality            
%                      constraints Ax <= b + Btheta.            
%                      Class: double                            
%          b           Right hand side of the inequality        
%                      constraints Ax <= b + Btheta.            
%                      Class: double                            
%          pB          Right hand side of the inequality        
%                      constraints for parameters Ax <= b +     
%                      Btheta.                                  
%                      Class: double                            
%          Ae          Linear part of the equality constraints  
%                      A_ex=b_e + Etheta .                      
%                      Class: double                            
%                      Default: []                              
%          be          Right hand side of the equality          
%                      constraints A_ex=b_e + Etheta .          
%                      Class: double                            
%                      Default: []                              
%          pE          Right hand side of the equality          
%                      constraints for parameters A_ex=b_e +    
%                      Etheta .                                 
%                      Class: double                            
%                      Default: []                              
%          lb          Lower bound for the decision variables x 
%                      >= lb.                                   
%                      Class: double                            
%                      Default: []                              
%          ub          Upper bound for the decision variables x 
%                      <= ub.                                   
%                      Class: double                            
%                      Default: []                              
%          Ath         Linear part of the inequality            
%                      constraints A_thetatheta <= b_theta.     
%                      Class: double                            
%                      Default: []                              
%          bth         Right hand side of the inequality        
%                      constraints A_thetatheta <= b_theta.     
%                      Class: double                            
%                      Default: []                              
%          vartype     A vector of strings determining the type 
%                      of the variable. Supported characters    
%                      are C (continuous), I (integer), B       
%                      (binary), N (semi-integer), S            
%                      (semi-continuous). Example: First        
%                      variable from three is binary, the rest  
%                      is continuous: vartype='BCC';            
%                      Class: char                              
%                      Default: "                               
%          solver      S string specifying which solver should  
%                      be called.                               
%                      Class: char                              
%                      Default: []                              
%          M           Linear matrix involved in LCP.           
%                      Class: double                            
%                      Default: []                              
%          q           Right hand side vector involved in LCP.  
%                      Class: double                            
%                      Default: []                              
%          Q           Linear matrix involved in parametric     
%                      formulation of LCP.                      
%                      Class: double                            
%                      Default: []                              
%          Matrices    Structure with the matrices defining     
%                      MPLP/MPQP problem as returned by         
%                      mpt_constructMatrices function. For      
%                      detailed description, see help           
%                      mpt_constructMatrices.                   
%                      Class: struct                            
%          constraints Yalmip set of constraints that formulate 
%                      the given problem.                       
%                      Class: lmi                               
%          objective   Yalmip variable that represents          
%                      objective value of the given problem.    
%                      Class: sdpvar                            
%          theta       Specification of parametric variables    
%                      involved in problem formulated in        
%                      Yalmip.                                  
%                      Class: sdpvar                            
%          x           Specification of decision variables      
%                      involved in problem formulated in        
%                      Yalmip.                                  
%                      Class: sdpvar                            
%                        
%  
%  
%  OUTPUT
%  ------
%     
%        
%          problem Object of the Opt class that defines the 
%                  optimization problem.                    
%                  Class: Opt                               
%                    
%  
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
%   mailto:colin.jones@epfl.ch 
%     
%    
%   (c) 2010-2013  Martin Herceg: ETH Zurich
%   mailto:herceg@control.ee.ethz.ch 
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
 
 
    properties (GetAccess=public, SetAccess=private)
        % LP/QP/pLP/pQP variables
        % J(th) = min 0.5*u'*H*u + (pF*th+f)'*u + th'*Y*th + C*th + c
        %         s.t.  A*u <= b  + pB*th
        %               Ae*u = be + pE*th
        %               lb  <= u <= ub
        %               Ath*th <= bth
        
        A  = []; b  = []; pB = [];
        Ae = []; be = []; pE = [];
        H  = []; f  = []; pF = [];
        lb = []; ub = [];
        Y = []; C = []; c = [];
        Ath = []; bth = [];
        
        % LCP variables
        % w - M*z = q + Q*th, w,z  >= 0, w'*z = 0
        M  = []; q  = []; Q  = [];
    end
    
    % public properties
    properties (GetAccess = public, SetAccess = public)
        Data  % any user-defined data
    end
    
    
    properties (SetAccess = private)
        n  = 0; % Problem dimension
        m  = 0; % Number of inequalities
        me = 0; % Number of equality constraints
        d  = 0; % Number of parameters
        
        % Problem types given as strings:
        % "LCP" - linear complementarity problem
        % "LP" - linear problem
        % "QP" - quadratic problem
        % "MILP" - mixed integer linear problem
        % "MIQP" - mixed integer quadratic problem

        problem_type = '';
        vartype = ''; % type of variables C-continuous, I-integer, B-binary, S-semicontinuous, N-semiinteger        
        isParametric = false;
        
        recover = []; % Mapping from solved problem to original problem data
        varOrder = [];
        Internal = [];
	end
	
	properties
		solver = ''
	end
    
	methods
		function set.solver(obj, new_solver)
			% Opt.solver setter
			
			if ~ischar(new_solver)
				error('The solver must be a string.');
			end
			obj.solver = upper(new_solver);
		end
	end
	
    methods(Access = public)
        
        % Constructor
        function opt = Opt(varargin)
            if nargin > 0
                if isa(varargin{1},'lmi') || isa(varargin{1},'constraint')
                    % convert from YALMIP data
                    opt = opt.setYalmipData(varargin{:});
                elseif isa(varargin{1},'struct')
                    % convert from MPT2.6 data
                    if isfield(varargin{1},'G') && isfield(varargin{1},'W') && isfield(varargin{1},'E')
                        opt = opt.setMPT26Data(varargin{1});
                    else
                        % call normal constructor
                        opt = opt.setData(varargin{:});
                    end
                else
                    % call normal constructor
                    opt = opt.setData(varargin{:});
                end
                
                % validate
                opt.validate;
                    
            else
                error('Empty problems not supported.');
            end
        end
        
        function K = feasibleSet(obj, arg)
            % Computes the feasible set of a given parametric problem
            %
            % For a parametric problem
            %    min  J(z, x)
            %    s.t. A*z <= b + pB*x
            %         Ae*z = be + pE*x
            % the feasible set K is the polyhedron
            %   K = { x | \exists z s.t. A*z<=b+pB*x, Ae*z=be+pE*x }
            %
            % This method implements two procedures to compute K:
            %   1) if K=prob.feasibleSet() is called, the feasible set is
            %      calculated by projection (can be expensive)
            %   2) if K=prob.feasibleSet(regions) is called with "regions"
            %      being the critical regions of the parametric solution,
            %      then K is constructed as follows:
            %         For each facet of each region do:
            %          a) compute the center of the facet
            %          b) take a small step accross the facet
            %          c) solve the problem for the new point
            %          d) if the problem is infeasible, add the facet to
            %             the feasible set 
            %
            % Syntax:
            %   K = prob.feasibleSet()
            %   K = prob.feasibleSet(method)
            %   K = prob.feasibleSet(regions)
            %
            % Inputs:
            %      prob: parametric problem as an Opt object
            %   regions: (optional) critical regions of the parametric
            %            solution
            %    method: (optional) string identificator of the projection
            %            method to use (see help Polyhedron/projection). By
            %            default we use the 'ifourier' method.
            %
            % Output:
            %         K: feasible set as a redundant H-polyhedron
            
            global MPTOPTIONS
            
            narginchk(1, 2);
            if ~obj.isParametric
                error('The problem is not parametric.');
            end
            if nargin==2
                if isa(arg, 'Polyhedron')
                    use_projection = false;
                    regions = arg;
                elseif isa(arg, 'char')
                    use_projection = true;
                    method = arg;
                else
                    error('The input must be either a string or an array of Polyhedron objects.');
                end
            else
                % default projection method
                use_projection = true;
                method = 'ifourier';
            end
            
            if use_projection
                % compute the feasible set via projection
                
                if isequal(lower(obj.problem_type), 'lcp') && ...
                        ~isfield(obj.Internal, 'constraints')
                    % This is a pLCP that was defined manually
                    %
                    % LCP constraints:
                    %   I*w - M*z = q + Q*x
                    %   w >= 0, z >= 0
                    %   Ath*x <= bth
                    
                    [nw, nz] = size(obj.M);
                    nb = length(obj.bth);
                    % y = [x; w; z]
                    Ae = [-obj.Q, eye(nw), -obj.M];
                    be = obj.q;
                    A = [obj.Ath, zeros(nb, nw+nz); ...
                        zeros(nw, obj.d), -eye(nw), zeros(nw, nz); ...
                        zeros(nz, obj.d+nw), -eye(nz)];
                    b = [obj.bth; zeros(nw+nz, 1)];
                    ZX = Polyhedron('Ae', Ae, 'be', be, 'A', A, 'b', b);

                else
                    % LP/QP constraints:
                    %     A*z <= b + pB*x
                    %    Ae*z == be + pE*x
                    %   Ath*x <= bth
                    
                    if isfield(obj.Internal, 'constraints')
                        % This is a pLCP that was generated by
                        % Opt/qp2lcp(). Convert it to the original pQP
                        % formulation to reduce the number of variables
                        obj = Opt(obj.Internal.constraints);
                    end
                    
                    if obj.me>0
                        % eliminate equalities, but do it on a copy of the object
                        obj = obj.copy();
                        obj.eliminateEquations();
                    end
                    
                    % construct the polyhedron 
                    % { (x, z) | A*x<=b+pB*x, Ath*x<=bth }
                    ZX = Polyhedron([-obj.pB, obj.A; ...
                        obj.Ath, zeros(length(obj.bth), obj.n)], ...
                        [obj.b; obj.bth]);
                end
                
                if ZX.isEmptySet()
                    % the feasible set is empty
                    K = Polyhedron.emptySet(obj.d);
                    return
                end
                
                % the feasible set is given as the projection of ZX onto
                % the parametric space
                K = ZX.projection(1:obj.d, method);
                
            else
                % construct the feasible set via critical regions
                
                % length of step over the facet
                step_size = MPTOPTIONS.rel_tol*10;
                Hf = [];
                t = tic;
                n_fails = 0;
                for i = 1:length(regions)
                    % for each region
                    if toc(t) > MPTOPTIONS.report_period
                        fprintf('progress: %d/%d\n', i, length(regions));
                        t = tic;
                    end
                    % make sure we have the minimal H-representation
                    % (we do redundancy elimination here since it can be
                    % costly in high dimensions; hence we want to give the
                    % uer an appropriate progress report)
                    regions(i).minHRep();
                    for j = 1:length(regions(i).b)
                        % for each facet of the i-th region:
                        % 1) compute a point on the facet
                        lpsol = regions(i).chebyCenter(j);
                        if lpsol.exitflag == MPTOPTIONS.OK
                            % 2) compute the point accross the j-th facet
                            x = lpsol.x + regions(i).A(j, :)'/norm(regions(i).A(j, :)')*step_size;
                            % 3) and solve the problem for the new point
                            qpsol = obj.solve(x);
                            if qpsol.exitflag ~= MPTOPTIONS.OK
                                % 4) infeasible => add this facet to the feasible set
                                Hf = [Hf; regions(i).H(j, :)];
                            end
                        else
                            % numerical problems
                            n_fails = n_fails + 1;
                        end
                    end
                end
                if n_fails > 0
                    fprintf('WARNING: failed to compute points on %d facet(s)\n', n_fails);
                end
                if isempty(Hf)
                    % numerical problems, return R^n
                    K = Polyhedron.fullSpace(regions(1).Dim);
                else
                    K = Polyhedron(Hf(:, 1:end-1), Hf(:, end));
                end
            end
        end
    end
    
%     methods
%         %% SET methods
%         % check if vartype is correct
%         function set.vartype(obj,val)
%             if ~isempty(val)
%                 if isnumeric(val)
%                     % convert to char if it is numeric
%                     val = char(val);
%                 end
%                 if ~isvector(val) || ~ischar(val)
%                     error('The argument must be a vector of strings.');
%                 end
%                 % checking if string is correct
%                 for i=1:length(val)
%                     if ~any(strcmpi(val(i),{'C','I','B','S','N'}))
%                         %C-continuous, I-integer, B-binary, S-semicontinuous, N-semiinteger
%                         error('Given string does not belong to gropu "C-continuous, I-integer, B-binary, S-semicontinuouos, N-semiinteger.');
%                     end
%                 end
%                 
%                 obj.vartype = val;
%             end            
%         end
%         
%     end
    
    methods (Access=protected)
        function U = copyElement(obj)
            % Creates a copy of the union
            %
            %   copy = U.copy()
            
            % Note: matlab.mixin.Copyable.copy() automatically properly
            % copies arrays and empty arrays, no need to do it here.
            % Moreover, it automatically creates the copy of proper class.
            U = copyElement@matlab.mixin.Copyable(obj);

            % we don't know what type of arguments can be put in the future
            % to Internal properties, so we check if a field contains a
            % "copy" method
            if isstruct(obj.Internal)
                nf = fieldnames(obj.Internal);
                for i=1:numel(nf)
                    x = obj.Internal.(nf{i});
                    if isobject(x) && ismethod(x, 'copy');
                        U.Internal.(nf{i}) = x.copy();
                    end
                end
            else
                if isobject(obj.Internal) && ismethod(obj.Internal,'copy');
                    U.Internal = obj.Internal.copy();
                end
            end
            % we don't know what type of arguments can be stored inside
            % Data, so we check if it contains a "copy" method - then use
            % it to create new object.
            if isstruct(obj.Data)
                nd = fieldnames(obj.Data);
                for i=1:numel(nd)
                    x = obj.Data.(nd{i});
                    if isobject(x) && ismethod(x,'copy');
                        U.Data.(nd{i}) = x.copy();
                    end
                end
            else
                if isobject(obj.Data) && ismethod(obj.Data,'copy');
                    U.Data = obj.Data.copy;
                end                
            end
            
        end
 
    end
end
