function toMatlab(obj, filename, function_to_export, varargin)
% Generates pure Matlab code for evaluation of a particular function
%
% Syntax:
% ---
%
%   obj.toMatlab(filename, function)
%
% Inputs:
% ---
%
%       obj: single BinTreePolyUnion object
%  filename: name of exported file (the '.m' extension will be added)
%  function: name of function to be exported
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
% Export the function 'primal' of polynion P to myprimal.m:
%   P.toMatlab('myprimal', 'primal')
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
narginchk(3, Inf);

%% validation
if sum([obj.Num])==0
	error('The object has no regions.');
elseif numel(obj)~=1
    error('Single BinTreePolyUnion object please.');
elseif ~all(obj.hasFunction(function_to_export))
	error('No such function "%s" in the object.', function_to_export);
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
fprintf(fid, '%%Evaluate function "%s" via a binary search tree\n', ...
	function_to_export);
fprintf(fid, '%% \n');
fprintf(fid, '%%  [value, region] = %s(x)\n', f_name);
fprintf(fid, '%%\n');
fprintf(fid, '%%See "help BinTreePolyUnion/toMatlab" for more information.\n');

% each splitting half-space is checked by H*[x; -1]<=0
fprintf(fid, 'x=x(:);xh=[x;-1];\n');

fprintf(fid, 'if numel(x)~=%d,error(''The input vector must have %d elements.'');end\n', ...
	obj(1).Dim, obj(1).Dim);

%% export data

% we capture the data of the array of unions as follows:
%   nz: range of the function to be exported
%   nx: number of parameters
%    T: matrix representing the binary tree (first nx+1 columns describe
%       the splitting hyperplanes, penultimate column is index of a node to
%       visit if a'*x<=b, the final column is pointer to a node if a'*x>0.
%       (negative index denotes leaf node)
%   fH: quadratic term of the function to be exported
%   fF: linear term of the function to be exported
%   fg: constant term of the function to be exported

nz = obj(1).Set(1).Functions(function_to_export).R;
data = struct('nz', nz, 'nx', obj(1).Dim, ...
	'T', mat2str(obj.Tree), ...
	'fH', [], 'fF', [], 'fg', []);

% prepare the data
nr = 0;
for i = 1:numel(obj)
	% extract parameters of the function to evaluate
	[fnz, fH, fF, fg] = get_function_data(obj(i), function_to_export);
	if any(fnz~=nz)
		error('In all regions the function "%s" must have range %d.', ...
			function_to_export, nz);
	end
	data.fH = [data.fH; fH];
	data.fF = [data.fF; fF];
	data.fg = [data.fg; fg];
end	
	
% write the data
fprintf(fid, 'T=%s;\n', data.T);
if ~isempty(data.fH)
	fprintf(fid, 'fH=%s;\n', mat2str(data.fH));
end
fprintf(fid, 'fF=%s;\n', mat2str(data.fF));
fprintf(fid, 'fg=%s;\n', mat2str(data.fg));

%% export code

% traverse the tree starting from the root node
fprintf(fid, 'i=1;nz=%d;nx=%d;z=NaN(%d,1);\n', data.nz, data.nx, data.nz);
fprintf(fid, 'while true\n');
fprintf(fid, 'h=T(i,1:%d)*xh;\n', data.nx+1);
fprintf(fid, 'if h<=0,i=T(i,%d);\n', data.nx+2);
fprintf(fid, 'else,i=T(i,%d);end\n', data.nx+3);
fprintf(fid, 'if i<0,i=-i;\n');
% region found, evaluate the function
fprintf(fid, 'z=fF((i-1)*nz+1:i*nz,:)*x+fg((i-1)*nz+1:i*nz);\n');
if ~isempty(data.fH)
    fprintf(fid, 'z=z+x''*fH((i-1)*nx+1:i*nx,:)*x;\n');
end
fprintf(fid, 'return\n');
% infeasible
fprintf(fid, 'elseif i==0,return;end\n');
fprintf(fid, 'end\n');

fclose(fid);

fprintf('Function "%s" evaluate via a binary tree exported to "%s".\n', ...
	function_to_export, full_path);

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
