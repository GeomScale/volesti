function lcp_compile
%
%  LCP_COMPILE: Compilation script for LCP solver 
%  ===============================================
%  
%  
%  SYNTAX
%  ------
%     
%      lcp_compile
%    
%  
%  DESCRIPTION
%  -----------
%     Assuming that a C-compiler is installed, this script compiles "lcp.c" and "lcp_sfun.c"
%  mex-functions to obtain matlab executables.
%  
%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010  Martin Herceg: ETH Zurich
%   mailto:herceg@control.ee.ethz.ch 
%  
%  
%  LICENSE
%  -------
%    
%    This program is free software; you can redistribute it and/or modify it under the terms of the GNU
%  General Public License as published by the Free Software Foundation; either version 2.1 of the
%  License, or (at your option) any later version.
%    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
%  even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%  General Public License for more details.
%    You should have received a copy of the GNU General Public License along with this library; if not,
%  write to the  Free Software Foundation, Inc.,  59 Temple Place, Suite 330,  Boston, MA 02111-1307
%  USA
 
 
 
 
 
 

%% compiles LCP solver as Matlab executable

% check if all files are present
fc = {'lcp.c', 'lcp_main.c', 'lcp_matrix.c', 'lumod_dense.c', 'lumod_dense.c',...
    'lcp_sfun.c', 'lcp.h', 'lcp_matrix.h', 'lumod_dense.h'};
fa = dir;
lf = {fa(:).name};
flags = false(size(fc));
for i=1:length(fc)
    if ~isempty(strmatch(fc{i},lf,'exact'))
        flags(i)=true;
    end
end
if ~all(flags)
    missing = fc(~flags);
    msg = 'lcp_compile: To compile LCP solver, following files "';
    for i=1:length(missing)
	    if i<length(missing)
	    	msg = strcat(msg,missing{i},'" "');
        else
	        msg = strcat(msg,missing{i});
        end
    end
    msg = strcat(msg,'" are missing in the current directory.');
    error(msg);
end

% compile
try
    disp('Compiling LCP solver ...');
    
    if isunix
        % compile under LINUX/MAC
        mex -largeArrayDims lcp.c lcp_main.c lcp_matrix.c lumod_dense.c  -lmwblas -lmwlapack
        mex -largeArrayDims lcp_sfun.c lcp_main.c lcp_matrix.c lumod_dense.c  -lmwblas -lmwlapack
    else
        % Windows compilation
        % check for a compiler
        cn = mex.getCompilerConfigurations;
        
        if ~isempty(regexpi(cn(1).Name,'microsoft'))
            dirstr='microsoft';
        elseif ~isempty(regexpi(cn(1).Name,'watcom'))
            dirstr='watcom';
        elseif ~isempty(regexpi(cn(1).Name,'lcc'))
            dirstr='lcc';
        else
            error('You don''t have a suitable C-compiler installed.');
        end

        str1 = 'mex -largeArrayDims lcp.c lcp_main.c lcp_matrix.c lumod_dense.c ';
        str2 = 'mex -largeArrayDims lcp_sfun.c lcp_main.c lcp_matrix.c lumod_dense.c ';
        
        if ~isempty(regexpi(cn(1).Name,'watcom'))
            lcplib = ['openwatcom',filesep,'lib',filesep,'lcp_mex.lib'];
            flib7s = ['openwatcom',filesep,'lib',filesep,'flib7s.lib'];
            eval([str1,' ',lcplib,' ',flib7s,' LINKFLAGS="$LINKFLAGS option nocaseexact"']);
            eval([str2,' ',lcplib,' ',flib7s,' LINKFLAGS="$LINKFLAGS option nocaseexact"']);
        else
            if strcmp(computer,'PCWIN')
                libs = ['''',matlabroot,'''',filesep,'extern',filesep,'lib',filesep,'win32',filesep,dirstr,filesep];
            elseif strcmp(computer,'PCWIN64')
                libs = ['''',matlabroot,'''',filesep,'extern',filesep,'lib',filesep,'win64',filesep,dirstr,filesep];
            end
            eval([str1,libs,'libmwblas.lib ',libs,'libmwlapack.lib']);
            eval([str2,libs,'libmwblas.lib ',libs,'libmwlapack.lib']);
        end
    end
    disp('Compilation of LCP solver finished.')
    disp('For help type "help lcp" and for a simple example type "lcp_test" at the matlab prompt.');
catch
    disp('Error appeared through compilation process.');
    disp('Check out you MEX settings, paths, and available libraries.');
end
