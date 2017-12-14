function toC(obj, function_name, file_name, tie_break_fcn)
%
%  TOC: Export of PWA/PWQ function to C-code 
%  ==========================================
%  
%  
%  SYNTAX
%  ------
%     
%      controller.toC('function')
%      controller.toC('function','filename')
%      controller.toC('function','filename','tie_break_fcn')
%    
%  
%  DESCRIPTION
%  -----------
%     The function toC() exports given piecewise affine (PWA) or piecewise
%  quadratic (PWQ) function to C-language including a sequential evaluation
%  routine. The PWA/PWQ function must be attached to the PolyUnion object.
%     If the file name is not provided, the default output name is mpt_getInput.
%    The export routine generates two files on the output: 
%    
%     - mpt_getInput.c - which contains the PWA/PWQ function including the
%     sequential search 
%     - mpt_getInput_mex.c - mex interface for evaluation in Matlab 
%    The file mpt_getInput_mex can be compiled inside Matlab and used for fast
%  evaluation of PWA/PWQ function. The compilation is invoked by mex routine as
%  follows:
%     mex mpt_getInput_mex
%    The PWA/PWQ function can be exported using the tie-break option if the
%  function is multiple valued. The tie-breaking option determines which value of
%  PWA/PWQ function will be evaluated based on the selecting the minimum in the
%  tie-breaking function. In this case, the tie-breaking function must be attached
%  to the PolyUnion object as well. If no tie-breaking function is provided, the
%  first found value in the sequential search of PWA/PWQ function is evaluated.
%    The function toC() can export the floating point numbers to single or double
%  precision. The default setting is double but this can be modified in global
%  options
%     modules.geometry.unions.PolyUnion.toC.
%  
%  INPUT
%  -----
%     
%        
%          function      Name of the attached PWA/PWQ function to 
%                        export.                                  
%                        Class: char                              
%          filename      Base name of the file to be generated.   
%                        Class: char                              
%          tie_break_fcn Name of the attached scalar PWA/PWQ      
%                        function to be used in tie-breaking      
%                        case.                                    
%                        Class: char                              
%                          
%  
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Martin Herceg: ETH Zurich
%   mailto:herceg@control.ee.ethz.ch 
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
 
 
global MPTOPTIONS

% argument checks
narginchk(2, 4);

if ~ischar(function_name)
    error('The function name must be given as a string.');
end
if nargin<3,
    file_name = 'mpt_getInput';
    tie_break_fcn='';
elseif nargin<4
    if isempty(file_name)
        file_name = 'mpt_getInput';
    end
    if ~ischar(file_name)
        error('The file name must be given as a string.');
    end
    tie_break_fcn='';
end
if ~ischar(tie_break_fcn)
    error('The function name must be given as a string.');
end

% get the short name if the full path is provided
[path_to_file,short_name] = fileparts(file_name);
if ~isempty(regexp(short_name,'[^a-zA-Z0-9_]', 'once'))
    error('The file name must contain only alphanumerical characters including underscore "_".');
end
% correct the file name with the proper extenstion ".c"
if isempty(path_to_file)
    file_name = [short_name,'.c'];
    mex_name = [short_name,'_mex.c'];
else
    file_name = [path_to_file, filesep, short_name,'.c'];
    mex_name = [path_to_file, filesep, short_name,'_mex.c'];
end

% check that all PolyUnions have the same list of functions
lF = obj(1).listFunctions;
for i=2:numel(obj)
    for j=1:numel(lF)
        if ~obj(i).hasFunction(lF{j})
            error('All objects in the array must have attached the same list of functions.');
        end
    end
end

% check that the PolyUnion is full-dimensional
for i=1:numel(obj)
   if ~obj(i).isFullDim       
       error('The export to C-code works only for full-dimensional partitions.');
   end    
end

% only PWA/PWQ functions are allowed
for i=1:numel(obj)
    % is the requested function present?
    if ~obj(i).hasFunction(function_name)
        error('No such function "%s" in the object.', function_name);
    end
    for j=1:obj(i).Num
        if ~(isa(obj(i).Set(j).Functions(function_name), 'AffFunction') || ...
            isa(obj(i).Set(j).Functions(function_name), 'QuadFunction') )
            error('Only quadratic and affine functions are supported.');
        end
        if ~isempty(tie_break_fcn)
            if ~(isa(obj(i).Set(j).Functions(tie_break_fcn), 'AffFunction') || ...
                    isa(obj(i).Set(j).Functions(tie_break_fcn), 'QuadFunction') )
                error('Only quadratic and affine tie-break functions are supported.');
            end
        end
    end
end


% check that all PolyUnions have requested functions and the range is the same
for i=1:numel(obj)
    R = obj(1).Set(1).Functions(function_name).R;
    for j=1:obj(i).Num
        if R~=obj(i).Set(j).Functions(function_name).R
            error('The requested function "%s" must have the same range in all objects in the array.',function_name);
        end
    end
    % is the tie-break function present?
    if ~isempty(tie_break_fcn)
        if ~obj(i).hasFunction(tie_break_fcn)
            error('No such function "%s" in the object.', tie_break_fcn);
        end
        for j=1:obj(i).Num
            if obj(i).Set(j).Functions(tie_break_fcn).R~=1
                error('The tie-break function "%s" must be scalar valued in all objects in the array.',tie_break_fcn);
            end
        end
    end
end


% only quadratic and affine functions are supported
isquadratic = isa(obj(1).Set(1).Functions(function_name), 'QuadFunction');
if ~isempty(tie_break_fcn)
    isquadraticTB = isa(obj(1).Set(1).Functions(tie_break_fcn), 'QuadFunction');
end

% single or double precision to export?
precision = MPTOPTIONS.modules.geometry.unions.PolyUnion.toC.precision;
precision = strtrim(lower(precision));
if isempty(precision)
    precision = 'double';
elseif ~isequal(precision,'single') && ~isequal(precision,'double')
    error('The specified precision in the option can be either "single" or "double".');
end
if isequal(precision,'single')
    precision = 'float';
end


% extract polyhedra with control law
Pn = [obj.Set];
nr = [obj.Num];
total_nr = sum(nr);

% extract hyperplane representation
An = cell(total_nr,1);
bn = cell(total_nr,1);
[An{:}]=deal(Pn.A);
[bn{:}]=deal(Pn.b);
if ~iscell(An),
    An = {An};
    bn = {bn};
end

% count number of constraints
nctotal = 0;
for ii=1:total_nr,
    nctotal = nctotal + size(Pn(ii).H,1);
end


% extract dimensions
nx = obj(1).Dim; % domain
nu = obj(1).Set(1).Functions(function_name).R; % range

% extract PWQ/PWA function
if isquadratic
    Hi = cell(total_nr,1);
end
Fi = cell(total_nr,1);
Gi = Fi;
for i=1:total_nr
    if isquadratic
        Hi{i}=Pn(i).Functions(function_name).H;
    end
    Fi{i}=Pn(i).Functions(function_name).F;
    Gi{i}=Pn(i).Functions(function_name).g;
end

% extract tie-break function
if ~isempty(tie_break_fcn)
    if isquadraticTB
        HTBi = cell(total_nr,1);
    end
    FTBi = cell(total_nr,1);
    GTBi = FTBi;
    for i=1:total_nr
        if isquadraticTB
            HTBi{i}=Pn(i).Functions(tie_break_fcn).H;
        end
        FTBi{i}=Pn(i).Functions(tie_break_fcn).F;
        GTBi{i}=Pn(i).Functions(tie_break_fcn).g;
    end
end

%% write mpt_getInput.c

fid = fopen(file_name, 'w');
if fid<0,
    error('Cannot open file "%s" for writing!',file_name);
end

fprintf(fid,'/* The function for evaluation of a piecewise affine control law associated\n');
fprintf(fid,'   to a given state X using sequential search.\n\n');
fprintf(fid,'  Usage:\n   unsigned long region = %s( %s *X, %s *U)\n\n',short_name, precision, precision);

header = {''    
'   where X is a vector of dimension MPT_DOMAIN and U is a vector of dimension '
'   MPT_RANGE for PWA functions ( 1 for PWQ functions). The output variable "region"'
'   indicates index of a region where the point X is located. If "region" index is'
'   smaller than 1 (region < 1), there is no control law associated to'
'   a given state.'
''
'   Please note that all code in this file is provided under the terms of the'
'   GNU General Public License, which implies that if you include it directly'
'   into your commercial application, you will need to comply with the license.'
'   If you feel this is not a good solution for you or your company, feel free '
'   to contact the author:'
'       michal.kvasnica@stuba.sk'
'    to re-license this specific piece of code to you free of charge.'
'*/'
''
'/* Copyright (C) 2005 by Michal Kvasnica (michal.kvasnica@stuba.sk) '
'   Revised in 2012-2013 by Martin Herceg (herceg@control.ee.ethz.ch)    '
'*/'
''
''
'/*  This program is free software; you can redistribute it and/or modify'
'    it under the terms of the GNU General Public License as published by'
'    the Free Software Foundation; either version 2 of the License, or'
'    (at your option) any later version.'
''
'    This program is distributed in the hope that it will be useful,'
'    but WITHOUT ANY WARRANTY; without even the implied warranty of'
'    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
'    GNU General Public License for more details.'
''
'    You should have received a copy of the GNU General Public License'
'    along with this program; if not, write to the Free Software'
'    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.'
'*/ '
''};

% write the header
for i=1:numel(header)
    fprintf(fid,[header{i},'\n']);
end

fprintf(fid,'/* Generated on %s by MPT %s */ \n\n',datestr(now), MPTOPTIONS.version);

if ~isempty(tie_break_fcn)
    fprintf(fid, '#include <float.h>\n\n');
end
fprintf(fid, '#define MPT_NR %d\n', total_nr);
fprintf(fid, '#define MPT_DOMAIN %d\n', nx);
if ~isquadratic
    fprintf(fid, '#define MPT_RANGE %d\n', nu);
end
fprintf(fid, '#define MPT_ABSTOL %e\n', MPTOPTIONS.abs_tol);


% write inequality constraints A*x <= b for each polytope
ctr = 0;
fprintf(fid, '\nstatic %s MPT_A[] = {\n', precision);
for ii = 1:total_nr,
    Ai = An{ii};
    nc = size(Ai, 1);
    for jj = 1:nc,
        a = Ai(jj, :);
        for kk = 1:length(a),
            ctr = ctr + 1;
            if ctr<nctotal*nx,
                if isequal(precision,'float')
                    fprintf(fid, '%.7e,\t', a(kk));
                else
                    fprintf(fid, '%.14e,\t', a(kk));
                end
            else
                if isequal(precision,'float')
                    fprintf(fid, '%.7e ', a(kk));
                else
                    fprintf(fid, '%.14e ', a(kk));
                end
            end
            if mod(ctr, 5)==0,
                fprintf(fid, '\n');
            end
        end
    end
end
fprintf(fid, '};\n\n');

ctr = 0;
fprintf(fid, 'static %s MPT_B[] = {\n',precision);
for ii = 1:total_nr,
    bi = bn{ii};
    nc = size(bi, 1);
    for jj = 1:nc,
        ctr = ctr + 1;
        if ctr<nctotal,
            if isequal(precision,'float')
                fprintf(fid, '%.7e,\t', bi(jj));
            else
                fprintf(fid, '%.14e,\t', bi(jj));
            end
        else
            if isequal(precision,'float')
                fprintf(fid, '%.7e ', bi(jj));
            else
                fprintf(fid, '%.14e ', bi(jj));
            end
        end
        if mod(ctr, 5)==0,
            fprintf(fid, '\n');
        end
    end
end
fprintf(fid, '};\n\n');

fprintf(fid, 'static int MPT_NC[] = {\n');
for ii = 1:total_nr,
    if ii < total_nr,
        fprintf(fid, '%d,\t', size(Pn(ii).H,1));
    else
        fprintf(fid, '%d ', size(Pn(ii).H,1));
    end
    if mod(ii, 5)==0,
        fprintf(fid, '\n');
    end
end
fprintf(fid, '};\n\n');


% write quadratic, linear and affine terms f(x) = x'*H*x + F*x + G
if isquadratic
    % quadratic term H
    sub_write_matrix(Hi, 'MPT_H', fid, precision);
end
% linear term F
sub_write_matrix(Fi, 'MPT_F', fid, precision);
% affine term G
sub_write_matrix(Gi, 'MPT_G', fid, precision);
if ~isempty(tie_break_fcn)
    if isquadraticTB
        % quadratic term H
        sub_write_matrix(HTBi, 'MPT_HTB', fid, precision);
    end
    % linear term F
    sub_write_matrix(FTBi, 'MPT_FTB', fid, precision);
    % affine term G
    sub_write_matrix(GTBi, 'MPT_GTB', fid, precision);
end
    
fprintf(fid, '\n\n/* main evaluation function using sequential search */\n');
fprintf(fid,'static unsigned long %s( %s *X, %s *U)\n', short_name, precision, precision);

if isempty(tie_break_fcn)
    if isquadratic
        % evaluation of a quadratic function in the first found region 
        core = {''
            '{'
            '    int ix, jx, ic, nc, isinside;'
            '    unsigned long ireg, abspos, region;'
            sprintf('    %s hx, sx;', precision)
            ''
            '    abspos = 0;'
            '    region = 0;'
            ''
            '    /* initialize U to zero*/'
            '    U[0] = 0;'
            ''
            '    for (ireg=0; ireg<MPT_NR; ireg++) {'
            ''
            '        isinside = 1;'
            '        nc = MPT_NC[ireg];'
            '        for (ic=0; ic<nc; ic++) {'
            '            hx = 0;'
            '            for (ix=0; ix<MPT_DOMAIN; ix++) {'
            '                hx += MPT_A[abspos*MPT_DOMAIN+ic*MPT_DOMAIN+ix]*X[ix];'
            '            }'
            '            if ((hx - MPT_B[abspos+ic]) > MPT_ABSTOL) {'
            '                /* constraint is violated, continue with next region */'
            '                isinside = 0;'
            '                break;'
            '            } '
            '        }'
            '        if (isinside==1) {'
            '            /* state belongs to this region, extract control law and exit */'
            '            region = ireg + 1;'
            '            for (ix=0; ix<MPT_DOMAIN; ix++) {'
            '                sx = 0;'
            '                for (jx=0; jx<MPT_DOMAIN; jx++) {'
            '                    sx += MPT_H[ireg*MPT_DOMAIN*MPT_DOMAIN + ix*MPT_DOMAIN + jx]*X[jx];'
            '                }'
            '                U[0] += sx*X[ix];'
            '            }'
            '            for (ix=0; ix<MPT_DOMAIN; ix++) {'
            '                U[0] += MPT_F[ireg*MPT_DOMAIN + ix]*X[ix];'
            '            }'
            '            U[0] += MPT_G[ireg];'
            '            return region;'
            '        }'
            '        abspos = abspos + MPT_NC[ireg];'
            '    }'
            '    return region;'
            '}'
            ''
            ''};        
    else
        % evaluation of an affine function in the first found region
        core = {''
            '{'
            '    int ix, iu, ic, nc, isinside;'
            '    unsigned long ireg, abspos, region;'
            sprintf('    %s hx;', precision)
            ''
            '    abspos = 0;'
            '    region = 0;'
            ''
            '    /* initialize U to zero*/'
            '    for (iu=0; iu<MPT_RANGE; iu++) {'
            '        U[iu] = 0;'
            '    }'
            ''
            '    for (ireg=0; ireg<MPT_NR; ireg++) {'
            ''
            '        isinside = 1;'
            '        nc = MPT_NC[ireg];'
            '        for (ic=0; ic<nc; ic++) {'
            '            hx = 0;'
            '            for (ix=0; ix<MPT_DOMAIN; ix++) {'
            '                hx = hx + MPT_A[abspos*MPT_DOMAIN+ic*MPT_DOMAIN+ix]*X[ix];'
            '            }'
            '            if ((hx - MPT_B[abspos+ic]) > MPT_ABSTOL) {'
            '                /* constraint is violated, continue with next region */'
            '                isinside = 0;'
            '                break;'
            '            } '
            '        }'
            '        if (isinside==1) {'
            '            /* state belongs to this region, extract control law and exit */'
            '            region = ireg + 1;'
            '            for (iu=0; iu<MPT_RANGE; iu++) {'
            '                for (ix=0; ix<MPT_DOMAIN; ix++) {'
            '                    U[iu] = U[iu] + MPT_F[ireg*MPT_DOMAIN*MPT_RANGE + iu*MPT_DOMAIN + ix]*X[ix];'
            '                }'
            '                U[iu] = U[iu] + MPT_G[ireg*MPT_RANGE + iu];'
            '            }'
            '            return region;'
            '        }'
            '        abspos = abspos + MPT_NC[ireg];'
            '    }'
            '    return region;'
            '}'
            ''
            ''};
    end
else
    % evaluation using tie-break function
    core = {''
            '{'
            '    int ix, jx, ic, nc, isinside;'
            '    unsigned long ireg, abspos, iregmin, region;'
            sprintf('    %s hx, sx, obj, objmin;', precision)
            ''
            '    abspos = 0;'
            '    region = 0;'
            '    iregmin = 0;'
            ''
            '   /* initialize values of the tie-break function */'
            '    obj = 0;'
            };
    if isequal(precision,'float')
        core = [core;   sprintf('    objmin = FLT_MAX;')];
    else
        core = [core;   sprintf('    objmin = DBL_MAX;')];
    end
    core = [core; {        
            ''
            ''
            '    for (ireg=0; ireg<MPT_NR; ireg++) {'
            ''
            '        isinside = 1;'
            '        nc = MPT_NC[ireg];'
            '        for (ic=0; ic<nc; ic++) {'
            '            hx = 0;'
            '            for (ix=0; ix<MPT_DOMAIN; ix++) {'
            '                hx += MPT_A[abspos*MPT_DOMAIN+ic*MPT_DOMAIN+ix]*X[ix];'
            '            }'
            '            if ((hx - MPT_B[abspos+ic]) > MPT_ABSTOL) {'
            '                /* constraint is violated, continue with next region */'
            '                isinside = 0;'
            '                break;'
            '            } '
            '        }'
            '        if (isinside==1) {'
            '            /* state belongs to this region, evaluate the tie-breaking function */'
            '            obj = 0;'
            }];
    if isquadraticTB
        core = [core; {
            '            for (ix=0; ix<MPT_DOMAIN; ix++) {'
            '                sx = 0;'
            '                for (jx=0; jx<MPT_DOMAIN; jx++) {'
            '                    sx += MPT_HTB[ireg*MPT_DOMAIN*MPT_DOMAIN + ix*MPT_DOMAIN + jx]*X[jx];'
            '                }'
            '                obj += sx*X[ix];'
            '            }'
            }];
    end
    core = [core; {
            '            for (ix=0; ix<MPT_DOMAIN; ix++) {'
            '                obj += MPT_FTB[ireg*MPT_DOMAIN + ix]*X[ix];'
            '            }'
            '            obj += MPT_GTB[ireg];'
            '            if (obj<objmin) {'
            '                objmin = obj;'
            '                region = ireg + 1;'
            '                iregmin = ireg;'
            '            }'
            '        }'
            '        abspos = abspos + MPT_NC[ireg];'
            '    }'
            }];
    if isquadratic
        core = [core; {
            '    U[0] = 0;'
            '    for (ix=0; ix<MPT_DOMAIN; ix++) {'
            '        sx = 0;'
            '        for (jx=0; jx<MPT_DOMAIN; jx++) {'
            '             sx += MPT_H[iregmin*MPT_DOMAIN*MPT_DOMAIN + ix*MPT_DOMAIN + jx]*X[jx];'
            '        }'
            '        U[0] += sx*X[ix];'
            '    }'
            '    for (ix=0; ix<MPT_DOMAIN; ix++) {'
            '        U[0] += MPT_F[iregmin*MPT_DOMAIN + ix]*X[ix];'
            '    }'
            '    U[0] += MPT_G[iregmin];'
            }];
    else
        core = [core; {
            '    for (ix=0; ix<MPT_RANGE; ix++) {'
            '        sx = 0;'
            '        for (jx=0; jx<MPT_DOMAIN; jx++) {'
            '            sx += MPT_F[iregmin*MPT_DOMAIN*MPT_RANGE + ix*MPT_DOMAIN + jx]*X[jx];'
            '        }'
            '        U[ix] = sx + MPT_G[iregmin*MPT_RANGE + ix];'
            '    }'
            }];
    end
    core = [core; {
            '    return region;'
            '}'
            ''
            ''}];
    
end

% write the core
for i=1:numel(core)
    fprintf(fid,[core{i},'\n']);
end

fclose(fid);
fprintf('Output written to "%s".\n', file_name);


%% write mpt_getInput_mex.c

fdn = fopen(mex_name, 'w');
if fdn<0,
    error('Cannot open file "%s" for writing!',mex_name);
end

% write the mex-file for evaluation under Matlab
cmex={''
    '/*'
    '  Autogenerated C-mex file for evalution of explicit controllers.'
    ''
    '  This file is to be compiled under Matlab using the syntax'
    ''
    sprintf('   mex -largeArrayDims %s.c',[short_name,'_mex'])
    ''
    '  Usage in Matlab:'
    ''
    sprintf('   u = %s(x)', [short_name,'_mex'])
    ''
    '  where x is the vector of double for which to evaluate PWA/PWQ function.'
    ''
    ''
    '  Copyright 2013 by Martin Herceg, Automatic Control Laboratory,'
    '  ETH Zurich, herceg@control.ee.ethz.ch'
    ''
    ''
    '  This program is free software; you can redistribute it and/or modify'
    '  it under the terms of the GNU General Public License as published by'
    '  the Free Software Foundation; either version 2 of the License, or'
    '  (at your option) any later version.'
    ''
    '  This program is distributed in the hope that it will be useful,'
    '  but WITHOUT ANY WARRANTY; without even the implied warranty of'
    '  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
    '  GNU General Public License for more details.'
    ''
    '  You should have received a copy of the GNU General Public License'
    '  along with this program; if not, write to the Free Software'
    '  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.'
    '*/ '
    ''
    ''
    sprintf('/* Generated on %s by MPT %s */ \n\n',datestr(now), MPTOPTIONS.version)
    ''
    ''
    '#include <math.h>'
    '#include "mex.h"'
    sprintf('#include "%s.c"',short_name)
    ''
    '#ifndef MPT_RANGE'
    '#define MPT_RANGE 1'
    '#endif'
    ''
    ''
    'void mexFunction( int nlhs, mxArray *plhs[], '
    '		  int nrhs, const mxArray*prhs[] )'
    ''
    '{ '
    '    int j;'
    '    double *xin, *uout;'
    '    char msg[70];'
    '    unsigned long region;'
    ''
    '    mwSize M,N,D; '
    };
if isequal(precision,'float')
    cmex = [cmex; '    float x[MPT_DOMAIN], u[MPT_RANGE];'];
end
cmex = [cmex; {
    ''
    '    /* Check for proper number of arguments */'
    ''
    '    if (nrhs != 1) { '
    '	    mexErrMsgTxt("One input arguments required."); '
    '    } else if (nlhs > 1) {'
    '	    mexErrMsgTxt("Too many output arguments."); '
    '    } '
    '    if (mxIsEmpty(prhs[0]))'
    '        mexErrMsgTxt("The argument must not be empty.");'
    '	if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) )'
    '        mexErrMsgTxt("The argument must be real and of type DOUBLE.");'
    '    if ( mxGetNumberOfDimensions(prhs[0])!=2 )'
    '		mexErrMsgTxt("The argument must be a vector"); '
    ''
    '    M = mxGetM(prhs[0]);'
    '    N = mxGetN(prhs[0]);'
    '    if ( (M!=1) && (N!=1) )'
    '		mexErrMsgTxt("The argument must be a vector.");'
    ''
    '    /* dimension */'
    '    D = mxGetNumberOfElements(prhs[0]);'
    '    if (D!=MPT_DOMAIN) {'
    '        sprintf(msg, "The size of the input argument must be %%d x 1.",MPT_DOMAIN);'
    '        mexErrMsgTxt(msg);'
    '    }'
    ''
    '    /* get and verify input data */'
    '    xin = mxGetPr(prhs[0]);   '
    '    for (j=0; j<D; j++) {'
    '        if ( mxIsNaN(xin[j]) || mxIsInf(xin[j]) )'
    '            mexErrMsgTxt("No ''NaN'' or ''Inf'' terms are allowed.");'
    '    }'
    ''
    '    /* create a matrix for the output */ '
    '    plhs[0] = mxCreateDoubleMatrix( MPT_RANGE, 1, mxREAL);    '
    '    uout = mxGetPr(plhs[0]);'
    ''
    }];
if isequal(precision,'float')    
    cmex = [cmex; {
    '    /* recast double as float */'
    '    for (j=0; j<MPT_DOMAIN; j++){'
    '        x[j] = (float)xin[j];'
    '    }'
    ''
    '    /* Do the evaluation in a subroutine */'
    sprintf('    region = %s(x, u);', short_name)
    ''
    }];
else
    cmex = [cmex; {        
    ''
    '    /* Do the evaluation in a subroutine */'
    sprintf('    region = %s(xin, uout);', short_name)
    ''    
    }];
end
if isequal(precision,'float')
    cmex = [cmex; {
    '    /* cast as double for output */'
    '    for (j=0; j<MPT_RANGE; j++){'
    '        uout[j] = (double)u[j];'
    '    }'
    ''}];
end
    cmex = [cmex; {
    ''
    '    return;'
    ''
    '}'
    }];

% write the core
for i=1:numel(cmex)
    fprintf(fdn,[cmex{i},'\n']);
end

fclose(fdn);
fprintf('C-mex function written to "%s".\n', mex_name);


end


function sub_write_matrix(matrices, name, fid, precision)
%
% writes a cell of matrices to a file with a given name
%
% syntax:
%         matrix - a cell array with matrix data
%         name - name of the matrix given as string
%         fid - file identificator
%         precision - either "float" or "double"

nr = numel(matrices);
[nu,nx] = size(matrices{1});


nctotalh = nu*nx*nr;
ctr = 0;
fprintf(fid, '\n');
fprintf(fid, 'static %s %s[] = {\n', precision, name);
for ii = 1:nr,
    M = matrices{ii};
    for jj = 1:nu,
        h = M(jj, :);
        for kk = 1:nx,
            ctr = ctr + 1;
            if ctr<nctotalh,
                if isequal(precision,'float')
                    fprintf(fid, '%.7e,\t', h(kk));
                else
                    fprintf(fid, '%.14e,\t', h(kk));
                end
            else
                if isequal(precision,'float')
                    fprintf(fid, '%.7e ', h(kk));
                else
                    fprintf(fid, '%.14e ', h(kk));
                end
            end
            if mod(ctr, 5)==0,
                fprintf(fid, '\n');
            end
        end
    end
end
fprintf(fid, '};\n\n');


end
