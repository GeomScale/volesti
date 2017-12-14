function toMatlab(obj, filename, function_to_export, tiebreak)
% Generates pure Matlab code for evaluation of a particular function
%
% Syntax:
% ---
%
%   obj.toMatlab(filename, function, tiebreak)
%
% Inputs:
% ---
%
%       obj: single PolyUnion or an array thereof
%  filename: name of exported file (the '.m' extension will be added)
%  function: name of function to be exported
%  tiebreak: name of function to use to resolve tiebreaks
%            (use 'first-region' to break the sequential search once
%            the first region containing a given point is found)
%
% Output:
% ---
%
%  This method generates "filename.m" that can be called as follows:
%
%    [zopt, region] = filename(x)
%
%  where X is the vector of parameters at which FUNCTION is to
%  be evaluated, ZOPT is the value of the function at X, and REGION
%  is the index of the region that contains X.
%
%  If no region contains the point X, then ZOPT=NaN and REGION=0.
%
% Example:
% ---
%
% Export the function 'primal' of polynion P to myprimal.m using the
% function 'obj' as a tie breaker:
%   P.toMatlab('myprimal', 'primal', 'obj')
%
% Evaluate the 'primal' function at a particular point:
%   [value, region] = myprimal(x)
%
% Convert the Matlab code to a mex function using the Matlab Coder:
%   coder myprimal -args {zeros(nx, 1)}
%   [value, region] = myprimal_mex(x)
% where NX is the number of elements of X.

global MPTOPTIONS

%% parsing
narginchk(4, 4);

%% validation
if sum([obj.Num])==0
	error('The object has no regions.');
elseif ~all(obj.hasFunction(function_to_export))
	error('No such function "%s" in the object.', function_to_export);
elseif isempty(tiebreak)
	error('The tiebreak function must be specified.');
elseif ~isequal(tiebreak, 'first-region') && ...
		~all(obj.hasFunction(tiebreak))
	error('No such function "%s" in the object.', tiebreak);
elseif any(diff([obj.Dim])~=0)
	error('All unions must be in the same dimension.');
elseif ~(isa(obj(1).Set(1).Functions(function_to_export), 'AffFunction') ||...
	isa(obj(1).Set(1).Functions(function_to_export), 'QuadFunction') )
	error('Only affine and quadratic functions can be exported.');
end

%% export headers
[f_path, f_name, f_ext] = fileparts(filename);
% replace invalid characters by an underscore
f_name = regexprep(f_name, '[^a-zA-Z0-9_]', '_');
if isempty(f_path)
	full_path = [f_name '.m'];
else
	full_path = [f_path filesep f_name '.m'];
end
fid = fopen(full_path, 'w');
fprintf(fid, 'function [z,i] = %s(x) %%#codegen\n', f_name);
fprintf(fid, '%%Evaluate function "%s" with tiebreak "%s"\n', ...
	function_to_export, tiebreak);
fprintf(fid, '%% \n');
fprintf(fid, '%%  [value, region] = %s(x)\n', f_name);
fprintf(fid, '%%\n');
fprintf(fid, '%%See "help PolyUnion/toMatlab" for more information.\n');

% H-representation is checked by H*[x; -1]<=0, He*[x; -1]==0
fprintf(fid, 'x=x(:);xh=[x;-1];\n');

fprintf(fid, 'if numel(x)~=%d,error(''The input vector must have %d elements.'');end\n', ...
	obj(1).Dim, obj(1).Dim);

%% export data

% we capture the data of the array of unions as follows:
%   nz: range of the function to be exported
%   nx: number of parameters
%    H: matrix of inequality H-representations of each region
%   He: matrix of equality H-representations of each region
%   ni: number of inequality constraints for each region
%   ne: number of equality constraints for each region
%   fH: quadratic term of the function to be exported
%   fF: linear term of the function to be exported
%   fg: constant term of the function to be exported
%   tH: quadratic term of the tiebreak function
%   tF: linear term of the tiebreak function
%   tg: constant term of the tiebreak function

nz = obj(1).Set(1).Functions(function_to_export).R;
data = struct('nz', nz, 'nx', obj(1).Dim, ...
	'H', [], 'He', [], 'ni', [], 'ne', [], ...
	'fH', [], 'fF', [], 'fg', [], ...
	'tH', [], 'tF', [], 'tg', []);

% prepare the data
nr = 0;
for i = 1:numel(obj)
	
	nr = nr + obj(i).Num;
	
	% extract H-representation of each region
	H = obj(i).Set.forEach(@(x) x.H, 'UniformOutput', false);
	He = obj(i).Set.forEach(@(x) x.He, 'UniformOutput', false);
	% number of inequality/equality constraints for each region
	ni = obj(i).Set.forEach(@(x) size(x.H, 1), 'UniformOutput', false);
	ne = obj(i).Set.forEach(@(x) size(x.He, 1), 'UniformOutput', false);

	% concatenate cells into one large matrix
	data.H = [data.H; cat(1, H{:})];
	data.ni = [data.ni; cat(1, ni{:})];
	data.He = [data.He; cat(1, He{:})];
	data.ne = [data.ne; cat(1, ne{:})];

	% extract parameters of the functions
	[fnz, fH, fF, fg] = get_function_data(obj(i), function_to_export);
	if any(fnz~=nz)
		error('In all regions the function "%s" must have range %d.', ...
			function_to_export, nz);
	end
	data.fH = [data.fH; fH];
	data.fF = [data.fF; fF];
	data.fg = [data.fg; fg];

	if ~isequal(tiebreak, 'first-region')
		[tnz, tH, tF, tg] = get_function_data(obj(i), tiebreak);
		if any(tnz~=1)
			error('The tie breaker must be a scalar-valued function.');
		end
		data.tH = [data.tH; tH];
		data.tF = [data.tF; tF];
		data.tg = [data.tg; tg];
	end
end	
	
% write the data
fprintf(fid, 'nx=%d;nz=%d;\n', data.nx, data.nz);
if ~isempty(data.H)
	fprintf(fid, 'H=%s;\n', mat2str(data.H));
	fprintf(fid, 'ni=%s;\n', mat2str(cumsum([1;data.ni])));
end
if ~isempty(data.He)
	fprintf(fid, 'He=%s;\n', mat2str(data.He));
	fprintf(fid, 'ne=%s;\n', mat2str(cumsum([1;data.ne])));
end
if ~isempty(data.fH)
	fprintf(fid, 'fH=%s;\n', mat2str(data.fH));
end
fprintf(fid, 'fF=%s;\n', mat2str(data.fF));
fprintf(fid, 'fg=%s;\n', mat2str(data.fg));
if ~isempty(data.tH)
	fprintf(fid, 'tH=%s;\n', mat2str(data.tH));
end
if ~isempty(data.tF)
	fprintf(fid, 'tF=%s;\n', mat2str(data.tF));
end
if ~isempty(data.tg)
	fprintf(fid, 'tg=%s;\n', mat2str(data.tg));
end

%% export code
if ~isequal(tiebreak, 'first-region')
	fprintf(fid, 'tb=[];\n');
end

% loop through each region
fprintf(fid, 'for i=1:%d,\n', nr);

% is "x" contained in the region?
if ~isempty(data.H)
	fprintf(fid, 'if all(H(ni(i):ni(i+1)-1,:)*xh<=%.g);\n', MPTOPTIONS.abs_tol);
end
if ~isempty(data.He)
	fprintf(fid, 'if all(abs(He(ne(i):ne(i+1)-1,:)*xh<=%.g));\n', MPTOPTIONS.abs_tol);
end

% what should happen when region contains x?
if isequal(tiebreak, 'first-region')
	% return the optimizer
	fprintf(fid, 'z=fF((i-1)*nz+1:i*nz,:)*x+fg((i-1)*nz+1:i*nz);\n');
	if ~isempty(data.fH)
		fprintf(fid, 'z=z+x''*fH((i-1)*nx+1:i*nx,:)*x;\n');
	end
	fprintf(fid, 'return\n');
else
	% more difficult case: record value of the tiebreak function
	% (note that we require the tiebreak function to be scalar-valued)
	fprintf(fid, 'tv=tF(i,:)*x+tg(i);\n');
	if ~isempty(data.tH)
		fprintf(fid, 'tv=tv+x''*tH((i-1)*nx+1:i*nx,:)*x;\n');
	end
	fprintf(fid, 'tb=[tb;i,tv];\n');
end	

if ~isempty(data.He)
	fprintf(fid, 'end\n');% end if He*x==0
end
if ~isempty(data.H)
	fprintf(fid, 'end\n');% end if H*x<=0
end
fprintf(fid, 'end\n'); % end i=1:nr

% resolve ties if necessary
if ~isequal(tiebreak, 'first-region')
	% in which partition/region is the tiebreak value minimal?
	fprintf(fid, 'if ~isempty(tb)\n');
	fprintf(fid, '[~,j]=min(tb(:,end));\n');
	fprintf(fid, 'i=tb(j,1);\n');
	% evaluate the function in this region
	fprintf(fid, 'z=fF((i-1)*nz+1:i*nz,:)*x+fg((i-1)*nz+1:i*nz);\n');
	if ~isempty(data.fH)
		fprintf(fid, 'z=z+x''*fH((i-1)*nx+1:i*nx,:)*x;\n');
	end
	fprintf(fid, 'return\n');
	fprintf(fid, 'end\n');
end

% return indication of infeasibility
fprintf(fid, 'i=0;z=NaN(%d,1);\n', data.nz);

% footer
fprintf(fid, 'end\n');

fclose(fid);

fprintf('Function "%s" with tiebreak "%s" exported to "%s".\n', ...
	function_to_export, tiebreak, full_path);

end

%-----------------------------------------------------------
function [nz, Hfm, Fm, gm] = get_function_data(obj, funname)

fun = obj.Set(1).Functions(funname);
if ~(isa(fun, 'AffFunction') || isa(fun, 'QuadFunction'))
	error('Only affine and quadratic functions can be exported.');
end

% range of the function in each region
nz = obj.Set.forEach(@(x) x.Functions(funname).R);
F = obj.Set.forEach(@(x) x.Functions(funname).F, 'UniformOutput', false);
g = obj.Set.forEach(@(x) x.Functions(funname).g, 'UniformOutput', false);
Fm = cat(1, F{:});
gm = cat(1, g{:});
if isa(fun, 'QuadFunction')
	Hf = obj.Set.forEach(@(x) x.Functions(funname).H, 'UniformOutput', false);
	Hfm = cat(1, Hf{:});
else
	Hfm = [];
end

end
