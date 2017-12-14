function toPython(obj, filename, function_to_export, tiebreak)
% Generates Python code for evaluation of a particular function
%
% Syntax:
% ---
%
%   obj.toPython(filename, function, tiebreak)
%
% Inputs:
% ---
%
%       obj: single PolyUnion or an array thereof
%  filename: name of exported file (the '.py' extension will be added)
%  function: name of function to be exported
%  tiebreak: name of function to use to resolve tiebreaks
%            (use 'first-region' to break the sequential search once
%            the first region containing a given point is found)
%
% Output:
% ---
%
%  This method generates "filename.py" that can be called as follows:
%
%    zopt = filename(x)
%
%  where X is the vector of parameters at which FUNCTION is to
%  be evaluated, ZOPT is the value of the function at X, and REGION
%  is the index of the region that contains X.
%
%  If no region contains the point X, then ZOPT=NaN.

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
	full_path = [f_name '.py'];
else
	full_path = [f_path filesep f_name '.py'];
end
fid = fopen(full_path, 'w');

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
% fprintf(fid, 'x=x(:);xh=[x;-1];\n');

fprintf(fid, 'from numpy import *\n');
fprintf(fid, 'import math\n');
fprintf(fid, 'def %s(x):\n',f_name);
fprintf(fid, '\txh=matrix(%s,dtype=object)\n',xh_python(data.nx));
fprintf(fid, '\tfor i in range(0,len(x)):\n');
fprintf(fid, '\t\txh[i]=x[i]\n');


fprintf(fid, '\tnx=%d;nz=%d;nu=%d;\n', data.nx, data.nz, nz);
if ~isempty(data.H)
	fprintf(fid, '\tH=matrix(%s);\n', mat2str_python(data.H));
	fprintf(fid, '\tni=matrix(%s);\n', mat2str_python(cumsum([1;data.ni])));
end
if ~isempty(data.He)
	fprintf(fid, '\tHe=matrix(%s);\n', mat2str_python(data.He));
	fprintf(fid, '\tne=matrix(%s);\n', mat2str_python(cumsum([1;data.ne])));
end
if ~isempty(data.fH)
	fprintf(fid, 'fH=matrix(%s);\n', mat2str_python(data.fH));
end
fprintf(fid, '\tfF=matrix(%s);\n', mat2str_python(data.fF));
fprintf(fid, '\tfg=matrix(%s);\n', mat2str_python(data.fg));
if ~isempty(data.tH)
	fprintf(fid, '\ttH=matrix(%s);\n', mat2str_python(data.tH));
end
if ~isempty(data.tF)
	fprintf(fid, '\ttF=matrix(%s);\n', mat2str_python(data.tF));
end
if ~isempty(data.tg)
	fprintf(fid, '\ttg=matrix(%s);\n', mat2str_python(data.tg));
end

%% export code
if ~isequal(tiebreak, 'first-region')
	fprintf(fid, '\ttb=array([]);\n');
end

% loop through each region
fprintf(fid, '\tfor i in range(0,%d):\n', nr);

% is "x" contained in the region?
if ~isempty(data.H)
    fprintf(fid, '\t\tif (H[ni[i]-1:ni[i+1]-1,:]*xh<=1e-08).all():\n', MPTOPTIONS.abs_tol);
end
if ~isempty(data.He)
    fprintf(fid, '\t\tif (abs(He[ne[i]-1:ne[i+1]-1,:]*xh<=1e-08)).all()\n', MPTOPTIONS.abs_tol);

end


% what should happen when region contains x?
if isequal(tiebreak, 'first-region')
	% return the optimizer
	fprintf(fid, '\t\t\tz=fF[i*nz:(i+1)*nz,:]*x+fg[i*nz:(i+1)*nz]\n');
	if ~isempty(data.fH)
		fprintf(fid, '\t\t\tz=z+x.T*fH[(i*nx:(i+1)*nx,:]*x;\n');
	end
	fprintf(fid, '\t\t\treturn z[0:nu]\n');
else
	% more difficult case: record value of the tiebreak function
	% (note that we require the tiebreak function to be scalar-valued)
    fprintf(fid, '\t\t\ttv=tF[i,:]*x+tg[i];\n');
    if ~isempty(data.tH)
		fprintf(fid, '\t\t\ttv=tv+x.T*tH[i*nx:(i+1)*nx,:]*x;\n');
    end
    fprintf(fid, '\t\t\ttb=insert(tb,shape(tb),i);\n');
    fprintf(fid, '\t\t\ttb=insert(tb,shape(tb),squeeze(asarray(tv)));\n');
% 	fprintf(fid, 'tb=[tb;i,tv];\n');
end	

% resolve ties if necessary
if ~isequal(tiebreak, 'first-region')
	% in which partition/region is the tiebreak value minimal?
    fprintf(fid, '\tnum_reg=shape(tb);num_reg=num_reg[0]/2;\n');
    fprintf(fid, '\ttb=tb.reshape(num_reg,2);\n');
    fprintf(fid, '\tif tb.sum()!=0:\n');
    fprintf(fid, '\t\ti=int(tb[tb[:,1].argmin(),0]);\n');
    fprintf(fid, '\t\tz=fF[i*nz:(i+1)*nz,:]*x+fg[i*nz:(i+1)*nz];\n');
    fprintf(fid, '\t\treturn z[0:nu]\n');
    fprintf(fid, '\telse:\n');
fprintf(fid, '\t\tz=matrix(%s)\n',toPython_nan(obj(1).Dim));
    
    fprintf(fid, '\t\treturn z\n');
end

% return indication of infeasibility
if isequal(tiebreak, 'first-region')
fprintf(fid,'\t\tif i==%d:\n',nr-1);
fprintf(fid,'\t\t\tz=matrix(%s)\n',toPython_nan(obj(1).Dim));
fprintf(fid,'\t\t\treturn z');
end

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
%-----------------------------------------------------------

function out = mat2str_python(matrix)
a = '[';
for i = 1:size(matrix,1)
    b = '[';
    for j = 1:size(matrix,2)
            if j<size(matrix,2)
                b = strcat(b,mat2str(matrix(i,j)),',');
            elseif (i==size(matrix,1))&(j==size(matrix,2))
               b = strcat(b,mat2str(matrix(i,j)),']');      
            else
                b = strcat(b,mat2str(matrix(i,j)),'],');
            end            
    end
    a = strcat(a,b);
end
out = strcat(a,']');
end

%-----------------------------------------------------------

function out = toPython_nan(nu)
a = '[';
for i=1:nu
    if i<nu
    a = strcat(a,'[float(''NaN'')],');    
    else
    a = strcat(a,'[float(''NaN'')]');        
    end
end
out = strcat(a,']');
end

%-----------------------------------------------------------
function a = xh_python(k)
a = '[';
for i=1:k
   a = strcat(a,'[0],');
end
a = strcat(a,'[-1]]');
end


