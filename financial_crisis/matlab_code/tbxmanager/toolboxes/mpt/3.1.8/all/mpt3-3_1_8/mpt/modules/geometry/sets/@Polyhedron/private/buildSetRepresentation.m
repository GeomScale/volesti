function optMat = buildSetRepresentation(obj)
%
% Compute matrices representing the set for the purposes of
% optimization
%

% Prefer H-representation. The primary purpose of the toolbox is parametric
% optimization where polyhedra are generated in H-rep. If both reps are
% available, this is usually a consequence of an inadvert H-to-V
% conversion, e.g. in plotting. Such a conversion can be unreliable.
if obj.hasHRep
    % Set containment:
    %   H(:,1:end-1)  * x <= H(:,end)
    %   He(:,1:end-1) * x == He(:,end)
    %   -inf <= x <= inf
    optMat.type = 'H';
    if ~isempty(obj.H_int)
        optMat.A  = obj.H_int(:,1:end-1);
        optMat.b  = obj.H_int(:,end);
    else
        optMat.A  = zeros(0,obj.Dim);
        optMat.b  = zeros(0,1);
    end
    if ~isempty(obj.He_int)
        optMat.Ae = obj.He_int(:,1:end-1);
        optMat.be = obj.He_int(:,end);
    else
        optMat.Ae = zeros(0,obj.Dim);
        optMat.be = zeros(0,1);
    end
    optMat.lb = -inf*ones(obj.Dim,1);
    optMat.ub =  inf*ones(obj.Dim,1);

elseif obj.hasVRep
    % Set containment:
	%   x = lam'*V + gam'*R
	%   lam>=0, gam>=0, 1'*lam=1

	optMat.type = 'V';
    nV = size(obj.V_int,1); 
	nR = size(obj.R_int,1);

	% lam>=0, gam>=0
    optMat.lb = [-inf*ones(obj.Dim,1);zeros(nV+nR,1)];
    optMat.ub = [ inf*ones(obj.Dim,1);ones(nV,1);inf*ones(nR,1)];

	% no inequality constraints
	optMat.A = zeros(0, obj.Dim+nV+nR);
	optMat.b = zeros(0, 1);

	% equality constraints:
	%   x = lam'*V + gam'*R
	%   1'*lam=1	
    %                [  x]   
    %  [-I  V' R']   [lam]   [0]
    %  [ 0' 1' 0'] * [gam] = [1]
	optMat.Ae = [-eye(obj.Dim) obj.V_int' obj.R_int'];
	optMat.be = zeros(obj.Dim, 1);
	if nV>0
		% only include 1'*lam==1 if we actually do have vertices
		optMat.Ae = [optMat.Ae; zeros(1,obj.Dim) ones(1,nV) zeros(1,nR)];
		optMat.be = [optMat.be; 1];
	end

else
	% no representation is available => empty set
	optMat.type = 'H';
	optMat.A = zeros(0, obj.Dim);
	optMat.b = zeros(0, 1);
	optMat.Ae = zeros(0, obj.Dim);
	optMat.be = zeros(0, 1);
	optMat.lb = ones(obj.Dim, 1); % empty set
	optMat.ub = zeros(obj.Dim, 1);
end

end
