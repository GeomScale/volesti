function I = integrate(P, function_name)
% Integrates a given function over a polytope
%
%   I = P.integrate(myfun)
%
% Inputs:
%       P: Polyhedron object or an array thereof
%          (each element must be bounded)
%   myfun: String name of the function to integrate
%          (the function must be either affine or quadratic)
%
% Outputs:
%       I: Value of the integral
%
% Literature: 
%  Lasserre, J.B., Avrachenkov, K.E.: The Multi-Dimensional Version of
%  $\int_a^b x^p dx$
%
% http://www-sop.inria.fr/members/Konstantin.Avratchenkov/pubs/simplex.ps

narginchk(1, 2);
if nargin<2
	function_name = '';
end
[function_name, msg] = P.validateFunctionName(function_name);
error(msg); % the error is only thrown if msg is not empty

% deal with arrays
if numel(P) > 1
	I = P.forEach(@(x) x.integrate(function_name));
	return
end

% we need the V-representation (doing the computation here significantly
% speeds up Polyhedron/isBounded)
P.minVRep();

% validation
if ~P.isBounded()
	error('The polyhedron must be bounded.');
end
fun = P.Functions(function_name);
if ~( isa(fun, 'AffFunction') || isa(fun, 'QuadFunction') )
	error('Function "%s" must be either affine or quadratic.', function_name);
end
if fun.R~=1
	error('Range of function "%s" must be 1.', function_name);
end

% triangulate
T = P.triangulate();
n = P.Dim;
I = 0;

% Per Remark 2.2 and eq. (2.9) of the Lasserre's paper, intergal of a
% quadratic homogeneous function x'Qx over a simplex Delta=convh(s1, ...,
% sn) is equal to volume(Delta)*trace(Q*Qt), where Qt is constructed below.
%
% Note: the formula in (2.9) only holds for homogeneous quadratic
% functions.
for it = 1:numel(T)
	% vertices of the simplex, must be column-wise
	S = T(it).V'; 
	% volume of the simplex
	vol = T(it).volume();

	if isa(fun, 'QuadFunction')
		% integral of the quadratic term
		Qt = zeros(n);
		for i = 1:n
			for j = 1:n
				qs = 0;
				for k = 1:n+1
					for l = k:n+1
						qs = qs + 0.5*(S(i, k)*S(j, l) + S(j, k)*S(i, l));
					end
				end
				Qt(i, j) = 2/((n+1)*(n+2))*qs;
			end
		end
		I = I + vol*trace(fun.H*Qt);
	end
	
	% Integral of a linear function over a simplex is obtained per formula
	% (2.3) of the Lasserre's paper
	I = I + vol/nchoosek(n+1,n)*fun.F*sum(S, 2);
	
	% Finally, integral of a constant is simply the scaled volume of the
	% simplex
	I = I + vol*fun.g;
end

end
