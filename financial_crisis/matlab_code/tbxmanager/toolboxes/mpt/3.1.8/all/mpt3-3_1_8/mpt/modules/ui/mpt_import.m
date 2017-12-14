function out = mpt_import(varargin)
%
%  MPT_IMPORT: Converts sysStruct and probStruct into an MPT3 prediction model 
%  ============================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      model = mpt_import(sysStruct, probStruct)
%    
%  
%  DESCRIPTION
%  -----------
%     model = mpt_import(sysStruct, probStruct) converts MPT2-styled control
%  problem, described by sysStruct and probStruct structures, into an MPT3
%  prediction model. The model can then be used for synthesis of MPC controllers
%  via controller = MPCController(model).
%    The following list shows which fields of sysStruct and probStruct can be
%  automatically converted and what are their respective equivalents in MPT3: 
%    
%     - sysStruct.xmin: model.x.min 
%     - sysStruct.xmax: model.x.max 
%     - sysStruct.umin: model.u.min 
%     - sysStruct.umax: model.u.max 
%     - sysStruct.ymin: model.y.min 
%     - sysStruct.ymax: model.y.max 
%     - probStruct.Q: model.x.penalty 
%     - probStruct.R: model.u.penalty 
%     - probStruct.Qy: model.y.penalty 
%     - probStruct.P_N: model.x.terminalPenalty 
%     - probStruct.Tset: model.x.terminalSet 
%     - probStruct.xref: model.x.reference 
%     - probStruct.uref: model.u.reference 
%     - probStruct.yref: model.y.reference 
%     - probStruct.tracking: currently only tracking=0 is supported 
%     - probStruct.Nc: model.u.block 
%     - probStruct.sxmax: model.x.softMin and model.x.softMax 
%     - probStruct.sumax: model.u.softMin and model.u.softMax 
%     - probStruct.symax: model.y.softMin and model.y.softMax 
%     - probStruct.subopt_lev: currently only subopt_lev=0 is supported. Use ctrl =
%     EMinTimeController(model) to obtain a minimum-time controller 
%  
%  
%  INPUT
%  -----
%     
%        
%          sysStruct  System structure                         
%                     Class: struct                            
%          probStruct Problem structure                        
%                     Class: struct                            
%                       
%  
%  
%  OUTPUT
%  ------
%     
%        
%          model Prediction model                         
%                  
%  
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
 
 
if nargin==1 && isstruct(varargin{1})
	% import from sysStruct
	out = importSysStruct(varargin{1});

elseif nargin==2 && isstruct(varargin{1}) && isstruct(varargin{2})
    % import sysStruct
    
	sysStruct = varargin{1};
	probStruct = varargin{2};
	
	[sysStruct, probStruct] = mpt_verifySysProb(sysStruct, probStruct);
	out = importSysStruct(sysStruct);
	
	if iscell(probStruct.Q) || iscell(probStruct.R)
		error('Time varying penalty matrices not yet supported.');
	end
	if isfield(probStruct, 'Qd')
		error('probStruct.Qd cannot be imported automatically.');
	end
	if isfield(probStruct, 'Qz')
		error('probStruct.Qz cannot be imported automatically.');
	end
	if isfield(probStruct, 'dref')
		error('probStruct.dref cannot be imported automatically.');
	end
	if isfield(probStruct, 'zref')
		error('probStruct.zref cannot be imported automatically.');
	end
	
	% add penalties
	out.u.penalty = create_penalty(probStruct.R, probStruct.norm);
	if isfield(probStruct, 'Qy')
		out.y.penalty = create_penalty(probStruct.Qy, probStruct.norm);
	else
		out.x.penalty = create_penalty(probStruct.Q, probStruct.norm);
	end
	
	% custom settings
	if ~iscell(sysStruct.A) && isfield(probStruct, 'Tconstraint') && ...
			probStruct.Tconstraint==1 && probStruct.norm==2 && ...
			(~isfield(probStruct, 'tracking') || probStruct.tracking==0)
		% add LQR terminal set/penalty for LTI systems with 2-norm cost
		fprintf('Computing LQR terminal set and penalty:\n');
		S = out.LQRSet();
		P = out.LQRPenalty();
		out.x.with('terminalPenalty');
		out.x.terminalPenalty = P;
		out.x.with('terminalSet');
		out.x.terminalSet = S;
	elseif ~isfield(probStruct, 'Qy') && probStruct.Tconstraint==1
		% add terminal state penalty to be consistent with MPT2
		out.x.with('terminalPenalty');
		out.x.terminalPenalty = create_penalty(probStruct.Q, probStruct.norm);
	end

	if isfield(probStruct, 'Nc')
		out.u.with('block');
		out.u.block.from = probStruct.Nc;
		out.u.block.to = probStruct.N;
	end
	if isfield(probStruct, 'P_N')
		if ~out.x.hasFilter('terminalPenalty')
			out.x.with('terminalPenalty');
		end
		out.x.terminalPenalty = create_penalty(probStruct.P_N, probStruct.norm);
	end
	if isfield(probStruct, 'Tset') && isfulldim(probStruct.Tset)
		out.x.with('terminalSet');
		out.x.terminalSet = Polyhedron('H', double(probStruct.Tset));
	end
	if isfield(probStruct, 'xref')
		out.x.with('reference');
		out.x.reference = probStruct.xref;
	end
	if isfield(probStruct, 'uref')
		out.u.with('reference');
		out.u.reference = probStruct.uref;
	end
	if isfield(probStruct, 'yref')
		out.y.with('reference');
		out.y.reference = probStruct.yref;
	end
	
	% soft state constraints
	if isfield(probStruct, 'Sx') && isfield(probStruct, 'sxmax')
		out.x.with('softMax');
		out.x.with('softMin');
		out.x.softMax.penalty = create_penalty(probStruct.Sx, probStruct.norm);
		out.x.softMin.penalty = create_penalty(probStruct.Sx, probStruct.norm);
		out.x.softMax.maximalViolation = probStruct.sxmax;
		out.x.softMin.maximalViolation = probStruct.sxmax;
	elseif isfield(probStruct, 'Sx')
		out.x.with('softMax');
		out.x.with('softMin');
		out.x.softMax.penalty = create_penalty(probStruct.Sx, probStruct.norm);
		out.x.softMin.penalty = create_penalty(probStruct.Sx, probStruct.norm);
	elseif isfield(probStruct, 'sxmax')
		out.x.with('softMax');
		out.x.with('softMin');
		out.x.softMax.maximalViolation = probStruct.sxmax;
		out.x.softMin.maximalViolation = probStruct.sxmax;
	end
	
	% soft input constraints
	if isfield(probStruct, 'Su') && isfield(probStruct, 'sumax')
		out.u.with('softMax');
		out.u.with('softMin');
		out.u.softMax.penalty = create_penalty(probStruct.Su, probStruct.norm);
		out.u.softMin.penalty = create_penalty(probStruct.Su, probStruct.norm);
		out.u.softMax.maximalViolation = probStruct.sumax;
		out.u.softMin.maximalViolation = probStruct.sumax;
	elseif isfield(probStruct, 'Su')
		out.u.with('softMax');
		out.u.with('softMin');
		out.u.softMax.penalty = create_penalty(probStruct.Su, probStruct.norm);
		out.u.softMin.penalty = create_penalty(probStruct.Su, probStruct.norm);
	elseif isfield(probStruct, 'sumax')
		out.u.with('softMax');
		out.u.with('softMin');
		out.u.softMax.maximalViolation = probStruct.sumax;
		out.u.softMin.maximalViolation = probStruct.sumax;
	end
	
	% soft output constraints
	if isfield(probStruct, 'Sy') && isfield(probStruct, 'symax')
		out.y.with('softMax');
		out.y.with('softMin');
		out.y.softMax.penalty = create_penalty(probStruct.Sy, probStruct.norm);
		out.y.softMin.penalty = create_penalty(probStruct.Sy, probStruct.norm);
		out.y.softMax.maximalViolation = probStruct.symax;
		out.y.softMin.maximalViolation = probStruct.symax;
	elseif isfield(probStruct, 'Sy')
		out.y.with('softMax');
		out.y.with('softMin');
		out.y.softMax.penalty = create_penalty(probStruct.Sy, probStruct.norm);
		out.y.softMin.penalty = create_penalty(probStruct.Sy, probStruct.norm);
	elseif isfield(probStruct, 'symax')
		out.y.with('softMax');
		out.y.with('softMin');
		out.y.softMax.maximalViolation = probStruct.symax;
		out.y.softMin.maximalViolation = probStruct.symax;
	end
	
	% tracking
	if isfield(probStruct, 'tracking') && probStruct.tracking>0

		if probStruct.tracking==1
			% deltau formulation
			out.u.with('deltaPenalty');
			if isfield(probStruct, 'Rdu')
				out.u.deltaPenalty = create_penalty(probStruct.Rdu, probStruct.norm);
			else
				out.u.deltaPenalty = create_penalty(probStruct.R, probStruct.norm);
			end
			out.u.penalty = [];
		end
		
		if isfield(probStruct, 'Qy')
			% output reference
			if ~out.y.hasFilter('reference')
				out.y.with('reference');
			end
			out.y.reference = 'free';
		else
			% state reference
			if ~out.x.hasFilter('reference')
				out.x.with('reference');
			end
			out.x.reference = 'free';
		end
	end

	
	if probStruct.subopt_lev==1
		obj.N = 1;
		fprintf('Use MinTimeController for minimum-time controllers.\n');
	elseif probStruct.subopt_lev==2
		error('M-step control strategy not yet supported.');
	end
	
    
elseif nargin==1 && isa(varargin{1}, 'mptctrl')
    % import mptctrl objects
    
	ctrl = varargin{1};
	sys = mpt_import(ctrl.details.origSysStruct, ctrl.details.origProbStruct);
	
	if isequal(ctrl.type, 'explicit')
		cs = struct(ctrl);
		cs.convex = true;
		data.system = sys;
		data.N = ctrl.probStruct.N;
		data.optimizer = mpt_mpsol2pu(cs);
		out = EMPCController(data);
	else
		out = MPCController(sys, ctrl.probStruct.N);
	end
	
elseif nargin==1 && isa(varargin{1}, 'polytope')
	% import polytope objects
	
	[H, K] = pelemfun(@double, varargin{1});
	out = [];
	for i = 1:length(H)
		out = [out, Polyhedron('A', H{i}, 'b', K{i})];
	end
	
else
	error('Unrecognized case.');
	
end
    
end

%---------------------------
function P = create_penalty(W, t)
% creates a weighted t-norm with weight W

switch t
	case 1,
		P = OneNormFunction(W);
	case Inf,
		P = InfNormFunction(W);
	case 2
		P = QuadFunction(W);
	otherwise
		error('Unrecognized norm type "%s".', num2str(t));
end

% P = Penalty(W, t);

end

%---------------------------
function out = importSysStruct(sysStruct)

sysStruct = mpt_verifySysStruct(sysStruct);

% first reject cases which we do not support yet
if isfield(sysStruct, 'noise') && isfulldim(sysStruct.noise)
	error('Systems with additive noise not yet supported.');
end
if isfield(sysStruct, 'Aunc')
	error('Systems with parametric uncertainties not yet supported.');
end

% import dynamics
if isfield(sysStruct, 'data') && isfield(sysStruct.data, 'MLD')
	% MLD system
	out = MLDSystem(sysStruct.data.MLD);
elseif iscell(sysStruct.A)
	% PWA system
	out = PWASystem(sysStruct);
else
	% LTI system
	out = LTISystem(sysStruct);
end

if isfield(sysStruct, 'Pbnd')
	out.x.with('initialSet');
	out.x.initialSet = toPolyhedron(sysStruct.Pbnd);
end

end
