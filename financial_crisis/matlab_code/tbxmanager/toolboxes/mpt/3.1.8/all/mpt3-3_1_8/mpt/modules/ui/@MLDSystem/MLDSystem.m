classdef MLDSystem < AbstractSystem
%
%  MLDSYSTEM: Represents mixed logical dynamical systems 
%  ======================================================
%  
%  
%  SYNTAX
%  ------
%     
%      sys = MLDSystem(source_file)
%    
%  
%  DESCRIPTION
%  -----------
%     This class represents mixed logical dynamical (MLD) systems of the form 
%                              +                                           
%                             x   = A x + B  u + B  d + B  z + B           
%                                          1      2      3      5          
%                                                                          
%                              y  = C x + D  u + D  d + D  z + D           
%                                          1      2      3      5          
%                                                                          
%                    E  d + E  z  <= E  u + E  x + E                       
%                     2      3        1      4      5                      
%     where x  is the state vector, u  is the vector of inputs, y  is the vector of
%  outputs, d  is the vector of auxiliary binary variables, and z  is the vector of
%  auxiliary real variables.
%    Each MLD system defines following properties: 
%    
%     - A, B1, B2, B3, B5: matrices of the state-update equation (read-only) 
%     - C, D1, D2, D3, D5: matrices of the output equation (read-only) 
%     - E1, E2, E3, E4, E5: matrices of the constraints (read-only) 
%     - nx, nu, ny, nd, nz: number of states, inputs, outputs, binaries, and aux
%     reals (read-only) 
%     - x: specifications of system's states (see help SystemSignal) 
%     - u: specifications of system's inputs (see help SystemSignal) 
%     - y: specifications of system's outputs (see help SystemSignal) 
%     - d: specifications of binary variables (see help SystemSignal) 
%     - d: specifications of auxiliary real variables (see help SystemSignal) 
%    MLD models are automatically created from a corresponding HYSDEL source file
%  upon calling the constructor as follows:
%     mldsys = MLDSystem('model.hys')
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2003-2013  Michal Kvasnica: STU Bratislava
%   mailto:michal.kvasnica@stuba.sk 
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
 
 
	properties(SetAccess=protected)
		% MLD model:
		% x^+ = A*x + B1*u + B2*d + B3*z + B5
		%   y = C*x + D1*u + D2*d + D3*z + D5
		% E2*d + E3*z <= E1*u + E4*x + E5
		
		A  % x^+ = A*x + B1*u + B2*d + B3*z + B5
		B1 % x^+ = A*x + B1*u + B2*d + B3*z + B5
		B2 % x^+ = A*x + B1*u + B2*d + B3*z + B5
		B3 % x^+ = A*x + B1*u + B2*d + B3*z + B5
		B5 % x^+ = A*x + B1*u + B2*d + B3*z + B5
		C %   y = C*x + D1*u + D2*d + D3*z + D5
		D1 %   y = C*x + D1*u + D2*d + D3*z + D5
		D2 %   y = C*x + D1*u + D2*d + D3*z + D5
		D3 %   y = C*x + D1*u + D2*d + D3*z + D5
		D5 %   y = C*x + D1*u + D2*d + D3*z + D5
		E1 % E2*d + E3*z <= E1*u + E4*x + E5
		E2 % E2*d + E3*z <= E1*u + E4*x + E5
		E3 % E2*d + E3*z <= E1*u + E4*x + E5
		E4 % E2*d + E3*z <= E1*u + E4*x + E5
		E5 % E2*d + E3*z <= E1*u + E4*x + E5
		S  % Full MLD structure returned from HYSDEL
		nz % Number of auxiliary real variables
		nd % Number of auxiliary binary variables
    end
    
    methods(Access=protected)
        
        function out = display_internal(obj)
            % Returns a string information about the system
            
            plural = @(s, n) [num2str(n) ' ' s repmat('s', 1, double(n~=1))];
            out = sprintf('%s with %s, %s, %s, %s, %s', class(obj), ...
                plural('state', obj.nx), ...
                plural('input', obj.nu), ...
                plural('output', obj.ny), ...
                plural('auxiliary real', obj.nz), ...
                plural('auxiliary bin', obj.nd));
        end
    end

    methods
        
        function obj = MLDSystem(source, parameters)
            % Constructor for MLD systems
            %
            % To import from a HYSDEL file:
            %    s = MLDSystem('hysdel_name.hys')
			%    s = MLDSystem('hysdel_name.hys', parameters)
            %
            % To import from an MLD structure:
            %    s = MLDSystem(S)

            if nargin == 0
                return
            end
            if nargin < 2
                parameters = [];
            end
            
            if isstruct(source)
                % import from MLD structure
                S = source;
                if isfield(S, 'MLDisvalid') && isfield(S, 'InputName')
                    % this is a HYSDEL3 structure, convert to HYSDEL2
                    % format
                    S = obj.h3toh2(S);
                end

            elseif ischar(source)
                % compile using HYSDEL
				S = obj.compile_hysdel(source, parameters);
            end

            % allow for MLD structures without affine terms
            if ~isfield(S, 'B5')
                S.B5 = zeros(S.nx, 1);
            end
            if ~isfield(S, 'D5')
                S.D5 = zeros(S.ny, 1);
            end
            
			obj.S = S;
			
			obj.A = S.A;
			obj.B1 = S.B1;
			obj.B2 = S.B2; 
			obj.B3 = S.B3; 
			obj.B5 = S.B5; 
			obj.C = S.C;
			obj.D1 = S.D1;
			obj.D2 = S.D2; 
			obj.D3 = S.D3; 
			obj.D5 = S.D5; 
			obj.E1 = S.E1;
			obj.E2 = S.E2; 
			obj.E3 = S.E3; 
			obj.E4 = S.E4;
			obj.E5 = S.E5;
			obj.nx = S.nx;
			obj.nu = S.nu;
			obj.ny = S.ny;
			obj.nz = S.nz;
			obj.nd = S.nd;

			% create variables
			x = SystemSignal(obj.nx);
            x.name = 'x';
            x.setKind('x');
            obj.addComponent('x', x);
			obj.x.min = S.xl;
			obj.x.max = S.xu;
			if S.nxb>0
				% last 'nxb' elements are binary
				obj.x.with('binary');
				obj.x.binary = S.nxr+1:S.nx;
			end
			
			u = SystemSignal(obj.nu);
            u.name = 'u';
            u.setKind('u');
            obj.addComponent('u', u);
			obj.u.min = S.ul;
			obj.u.max = S.uu;
			if S.nub>0
				% last 'nub' elements are binary
				obj.u.with('binary');
				obj.u.binary = S.nur+1:S.nu;
			end
			
			y = SystemSignal(obj.ny);
            y.name = 'y';
            y.setKind('y');
            obj.addComponent('y', y);
			obj.y.min = S.yl;
			obj.y.max = S.yu;
			if S.nyb>0
				% last 'nyb' elements are binary
				obj.y.with('binary');
				obj.y.binary = S.nyr+1:S.ny;
			end
			
			d = SystemSignal(obj.nd);
            d.name = 'd';
            d.setKind('d');
            obj.addComponent('d', d);
			obj.d.min = S.dl;
			obj.d.max = S.du;
            if obj.nd>0
                obj.d.with('binary');
                obj.d.binary = true;
            end
			
			z = SystemSignal(obj.nz);
            z.name = 'z';
            z.setKind('z');
            obj.addComponent('z', z);
			obj.z.min = S.zl;
			obj.z.max = S.zu;

        end
        
        function out = toPWA(obj)
            % Converts the MLD form to a PWA system
            
            sysStruct = mpt_pwa2sys(hys2pwa(obj.S), obj.S);
			
			% TODO: check that we correctly propagate binarity of states,
			% inputs and outputs
            out = PWASystem(sysStruct);
        end

        function C = constraints(obj)
            % convert MLD model at time step "k" into YALMIP constraints
            
            % constraints on variables
            C = constraints@AbstractSystem(obj);

            % get MLD variables
            x = obj.x.var;
            u = obj.u.var;
            d = obj.d.var;
            z = obj.z.var;
            y = obj.y.var;
            S = obj.S;

            % add the MLD constraints
            for k = 1:obj.Internal.system.N
                if ~isempty(x)
                    C = C + [ x(:, k+1) == S.A*x(:, k) + S.B1*u(:, k) + ...
                        S.B2*d(:, k) + S.B3*z(:, k) + S.B5 ];
                end
                if ~isempty(y)
                    C = C + [ y(:, k) == S.C*x(:, k) + S.D1*u(:, k) + ...
                        S.D2*d(:, k) + S.D3*z(:, k) + S.D5 ];
                end
                C = C + [ S.E2*d(:, k) + S.E3*z(:, k) <= ...
                    S.E1*u(:, k) + S.E4*x(:, k) + S.E5 ];
            end
		end
		
	end
	
	methods(Hidden)
		% implementation of abstract methods

		% no validation in these functions! it was already performed in
		% AbstractSystem/update() and output()

		function varargout = update_equation(obj, x0, u)
            % Evaluates the state-update and output equations and updates
            % the internal state of the system
            %
			% [xn, y, z, d] = obj.update_equation(x, u)
			
			global MPTOPTIONS
			if isempty(MPTOPTIONS)
				MPTOPTIONS = mptopt;
			end
			
			feasible = false;
			if size(obj.E1, 1)==0
				% no constraints, compute 'xn' and 'y' directly
				xn = obj.A*x0 + obj.B1*u + obj.B5;
				y = obj.C*x0 + obj.D1*u + obj.D5;
				d = zeros(obj.nd, 1);
				z = zeros(obj.nz, 1);
				feasible = true;
				
			elseif obj.nz==0 && obj.nd==0
				% we have constraints, but not 'd' and 'z' variables, just
				% check constraints and compute 'xn' and 'y' directly
				b = obj.E1*u + obj.E4*x0 + obj.E5;
				d = [];
				z = [];
				if any(b < -MPTOPTIONS.zero_tol)
                    feasible = false;
					warning('MLD constraints lead to infeasible or unbounded solution.');
				else
					xn = obj.A*x0 + obj.B1*u + obj.B5;
					y = obj.C*x0 + obj.D1*u + obj.D5;
					feasible = true;
				end
				
			else
				% constraints + 'd' or 'z' variables, solve a feasibility
				% MILP
				
				% just feasibility
				cost = zeros(obj.nd + obj.nz, 1);
				
				% inequality constraints
				A = [obj.E2, obj.E3];
				b = obj.E1*u + obj.E4*x0 + obj.E5;
				
				% indicate binary variables (deltas are first)
				vartype = repmat('C', 1, obj.nd+obj.nz);
				vartype(1:obj.nd) = 'B';
				
				% lower and upper bounds of "d" and "z" variables
				lb = [zeros(obj.nd, 1); obj.S.zl];
				ub = [ones(obj.nd, 1); obj.S.zu];
				
				% formulate and solve the feasibility problem
				xopt = Opt('A', A, 'b', b, 'f', cost, 'vartype', vartype).solve();
				feasible = xopt.exitflag==MPTOPTIONS.OK;
				d = xopt.xopt(1:obj.nd);
				z = xopt.xopt(obj.nd+1:end);
				if feasible
                    xn = obj.A*x0 + obj.B5;
                    if ~isempty(obj.B1);
                        xn = xn + obj.B1*u;
                    end
                    if ~isempty(obj.B2)
                        xn = xn + obj.B2*d;
                    end
                    if ~isempty(obj.B3)
                        xn = xn + obj.B3*z;
                    end
                    y = obj.C*x0 + obj.D5;
					if ~isempty(obj.D1)
                        y = y + obj.D1*u;
                    end
                    if ~isempty(obj.D2)
                        y = y + obj.D2*d;
                    end
                    if ~isempty(obj.D3)
                        y = y + obj.D3*z;
                    end
				end
			end
			
			if ~feasible
				xn = NaN(obj.nx, 1);
				y = NaN(obj.ny, 1);
				d = NaN(obj.nd, 1);
				z = NaN(obj.nz, 1);
			end
			varargout{1} = xn;
			varargout{2} = y;
			varargout{3} = z;
			varargout{4} = d;
			
		end
        
		function y = output_equation(obj, x, u)
            % Evaluates the output equation
        
			warning('Computing the output for MLD systems is costly, use the "update()" method instead.');
			[~, y] = obj.update_equation(x, u);
		end
		
		function out = has_feedthrough(obj)
			% feedthrough indication. must return true if the system has
			% direct feedthrough, false otherwise

			out = true;
		end

	end
	
	
	methods(Static)
        
        function S = sort_symtable(S)
            % sort variables in symtable according to their index
            
            index = zeros(1, length(S.symtable));
            for i = 1:length(S.symtable)
                if ~isfield(S.symtable{i}, 'index')
                    index(i) = -1;
                elseif S.symtable{i}.type == 'b'
                    index(i) = 10000+S.symtable{i}.index;
                else
                    index(i) = S.symtable{i}.index;
                end
            end
            [a, b] = sort(index);
            S.symtable = S.symtable(b);
        end
        
        function [S, hysname] = compile_hysdel(hysname, params)
            % compile a given hysdel file
            
            % get the true hysdel model name
            [p, short_hysname, e] = fileparts(hysname);
            if ~isempty(p)
                error('Hysdel model "%s" must be in the current working directory', hysname);
            end
            h_name = [short_hysname '.hys'];
            m_name = [short_hysname '.m'];

            % first remove the autogenerated m-file
            if exist(['.' filesep m_name], 'file')
                delete(m_name);
            end
            
            % compile the hysdel source
            fprintf('Compiling file "%s"...\n', h_name);
            T = evalc('hysdel(h_name, '''', ''-5'');');
            if ~isempty(strfind(T, 'Error:'))
                disp(T);
                error('HYSDEL error.');
            end
            
            % get the MLD structure by executing the m-file
            pause(0.1);
            rehash;
            if ~exist(m_name, 'file')
                disp(T);
                error('HYSDEL error.');
            end
            run(short_hysname);
            delete(m_name);
            
            % sort the symtable
            S = MLDSystem.sort_symtable(S);
        end
        
        function MLDold = h3toh2(MLD)
            % HYSDEL3-to-HYSDEL2 conversion of the MLD structure
            
            % reorder variables such that real components are in frton
            % of binaries
            xo = [MLD.j.xr(:); MLD.j.xb];
            uo = [MLD.j.ur(:); MLD.j.ub];
            yo = [MLD.j.yr(:); MLD.j.yb];
            
            % state-update equation
            MLDold.A = MLD.A(xo, xo);
            MLDold.B1 = MLD.Bu(xo, uo);
            MLDold.B2 = MLD.Baux(xo, MLD.j.d);
            MLDold.B3 = MLD.Baux(xo, MLD.j.z);
            MLDold.B5 = MLD.Baff(xo);
            
            % output equation
            MLDold.C = MLD.C(yo, xo);
            MLDold.D1 = MLD.Du(yo, uo);
            MLDold.D2 = MLD.Daux(yo, MLD.j.d);
            MLDold.D3 = MLD.Daux(yo, MLD.j.z);
            MLDold.D5 = MLD.Daff(yo);
            
            % inequalities
            MLDold.E1 = -MLD.Eu(MLD.j.ineq, uo);
            MLDold.E2 = MLD.Eaux(MLD.j.ineq, MLD.j.d);
            MLDold.E3 = MLD.Eaux(MLD.j.ineq, MLD.j.z);
            MLDold.E4 = -MLD.Ex(MLD.j.ineq, xo);
            MLDold.E5 = MLD.Eaff(MLD.j.ineq);
            
            % add equalities as double-sided inequalities
            MLDold.E1 = [MLDold.E1; ...
                -MLD.Eu(MLD.j.eq, uo); MLD.Eu(MLD.j.eq, uo)];
            MLDold.E2 = [MLDold.E2; ...
                MLD.Eaux(MLD.j.eq, MLD.j.d); -MLD.Eaux(MLD.j.eq, MLD.j.d)];
            MLDold.E3 = [MLDold.E3; ...
                MLD.Eaux(MLD.j.eq, MLD.j.z); -MLD.Eaux(MLD.j.eq, MLD.j.z)];
            MLDold.E4 = [MLDold.E4; ...
                -MLD.Ex(MLD.j.eq, xo); MLD.Ex(MLD.j.eq, xo)];
            MLDold.E5 = [MLDold.E5; ...
                MLD.Eaff(MLD.j.eq); -MLD.Eaff(MLD.j.eq)];
            
            % dimensions
            MLDold.nx = MLD.nx;
            MLDold.nxr = length(MLD.j.xr);
            MLDold.nxb = length(MLD.j.xb);
            MLDold.nu = MLD.nu;
            MLDold.nur = length(MLD.j.ur);
            MLDold.nub = length(MLD.j.ub);
            MLDold.ny = MLD.ny;
            MLDold.nyr = length(MLD.j.yr);
            MLDold.nyb = length(MLD.j.yb);
            MLDold.nd = length(MLD.j.d);
            MLDold.nz = length(MLD.j.z);
            MLDold.ne = size(MLDold.E5, 1);
            
            % upper/lower bounds
            MLDold.xl = MLD.xl(xo);
            MLDold.xu = MLD.xu(xo);
            MLDold.ul = MLD.ul(uo);
            MLDold.uu = MLD.uu(uo);
            MLDold.yl = MLD.yl(yo);
            MLDold.yu = MLD.yu(yo);
            MLDold.dl = MLD.wl(MLD.j.d);
            MLDold.du = MLD.wu(MLD.j.d);
            MLDold.zl = MLD.wl(MLD.j.z);
            MLDold.zu = MLD.wu(MLD.j.z);
            MLDold.symtable = [];
            MLDold.MLDisvalid = 1;
        end
    end

end
