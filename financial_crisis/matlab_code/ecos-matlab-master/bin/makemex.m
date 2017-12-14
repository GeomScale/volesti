function makemex(what)
% Matlab MEX makefile for ECOS.
%
%    MAKEMEX(WHAT) is a make file for ECOS, the embedded conic solver. It
%    builds ECOS and its components from source for the current platform
%    for which this script is invoked. WHAT is a cell array of strings,
%    with the following options (zero, one or more are possible):
%
%    {}, '' (empty string) or 'all': build all components and link.
%
%    'ldl': build the sparse LDL solver
%
%    'amd': builds the AMD package
%
%    'ecos': builds the core ECOS solver
%
%    'ecos_bb': builds the branch and bound module for mixed-binary
%    problems
%
%    'ecosmex': builds the ECOS mex interface - this involves linking of
%    the packages LDL, AMD, ECOS and ECOS_BB so these must have been built 
%    before a "MAKEMEX ecosmex" will terminate successfully.
%
%
%    In addition, you can invoke MAKEMEX as follows:
%
%        makemex clean - delete all object files (.o and .obj)
%        makemex purge - same as above, and also delete the mex files for
%                        the specific platform the script is executed on. 
%
% (c) A. Domahidi, ETH Zurich & embotech GmbH, Zurich, Switzerland, 2012-15.  


if( nargin == 0 )
    what = {'all'};
elseif( isempty(strfind(what, 'all'))     && ...
        isempty(strfind(what, 'ldl'))     && ...
        isempty(strfind(what, 'amd'))     && ...
        isempty(strfind(what, 'clean'))   && ...
        isempty(strfind(what, 'purge'))   && ...
        isempty(strfind(what, 'ecos'))    && ...
        isempty(strfind(what, 'ecos_bb')) && ...
        isempty(strfind(what, 'ecosmex')) )
    fprintf('No rule to make target "%s", exiting.\n', what);
end


if (~isempty (strfind (computer, '64')))
    d = '-largeArrayDims -DDLONG -DLDL_LONG';    
else
    d = '-DDLONG -DLDL_LONG';    
end


%% ldl solver
if( any(strcmpi(what,'ldl')) || any(strcmpi(what,'all')) )
    fprintf('Compiling LDL...');
    cmd = sprintf('mex -c -O %s -I../ecos/external/ldl/include -I../ecos/external/SuiteSparse_config ../ecos/external/ldl/src/ldl.c', d);
    eval(cmd);
    fprintf('\t\t\t[done]\n');
end


%% amd (approximate minimum degree ordering)
if( any(strcmpi(what,'amd')) || any(strcmpi(what,'all')) )
    fprintf('Compiling AMD...');
    i = sprintf ('-I../ecos/external/amd/include -I../ecos/external/SuiteSparse_config') ;
    cmd = sprintf ('mex -c -O %s -DMATLAB_MEX_FILE %s', d, i) ;
    files = {'amd_order', 'amd_dump', 'amd_postorder', 'amd_post_tree', ...
             'amd_aat', 'amd_2', 'amd_1', 'amd_defaults', 'amd_control', ...
             'amd_info', 'amd_valid', 'amd_global', 'amd_preprocess' } ;
    for i = 1 : length (files)
        cmd = sprintf ('%s ../ecos/external/amd/src/%s.c', cmd, files {i}) ;
    end
    eval(cmd);
    fprintf('\t\t\t[done]\n');
end


%% ecos (core IPM solver)
if( any(strcmpi(what,'ecos')) || any(strcmpi(what,'all')) ) 
    fprintf('Compiling ecos...');
    i = sprintf('-I../ecos/include -I../ecos/external/SuiteSparse_config -I../ecos/external/ldl/include -I../ecos/external/amd/include');
    cmd = sprintf('mex -c -O %s -DMATLAB_MEX_FILE -DCTRLC=1 %s', d, i);
    files = {'ecos', 'kkt', 'cone', 'spla', 'ctrlc', 'timer', 'preproc', 'splamm', 'equil', 'wright_omega', 'expcone'};
    for i = 1 : length (files)
        cmd = sprintf ('%s ../ecos/src/%s.c', cmd, files {i}) ;
    end
    eval(cmd);    
    fprintf('\t\t\t[done]\n');
end


%% ecos_bb (branch and bound module)
if( any(strcmpi(what,'ecos_bb')) || any(strcmpi(what,'all')) ) 
    fprintf('Compiling ecos_bb...');
    i = sprintf('-I../ecos/include -I../ecos/external/SuiteSparse_config -I../ecos/external/ldl/include -I../ecos/external/amd/include');
    cmd = sprintf('mex -c -O %s -DMATLAB_MEX_FILE -DCTRLC=1 %s', d, i);
    files_bb = {'ecos_bb', 'ecos_bb_preproc'};
    for i = 1 : length (files_bb)
        cmd = sprintf ('%s ../ecos/ecos_bb/%s.c', cmd, files_bb {i}) ;
    end
    eval(cmd);    
    fprintf('\t\t\t[done]\n');
end


%% ecos_mex
if( any(strcmpi(what,'ecosmex')) || any(strcmpi(what,'all')) )
    fprintf('Compiling ecos_mex...');
    cmd = sprintf('mex -c -O -DMATLAB_MEX_FILE -DCTRLC=1 -DMEXARGMUENTCHECKS %s -I../ecos/include -I../ecos/external/SuiteSparse_config -I../ecos/external/ldl/include -I../ecos/external/amd/include ../src/ecos_mex.c', d);
    eval(cmd);
    fprintf('\t\t\t[done]\n');
    fprintf('Linking...     ');
    clear ecos
    if( ispc )
        cmd = sprintf('mex %s -lut amd_1.obj amd_2.obj amd_aat.obj amd_control.obj amd_defaults.obj amd_dump.obj amd_global.obj amd_info.obj amd_order.obj amd_post_tree.obj amd_postorder.obj amd_preprocess.obj amd_valid.obj ldl.obj kkt.obj preproc.obj spla.obj cone.obj ecos.obj ctrlc.obj timer.obj splamm.obj equil.obj ecos_mex.obj ecos_bb.obj ecos_bb_preproc.obj wright_omega.obj expcone.obj -output "ecos"', d);
        eval(cmd);    
    elseif( ismac )
        cmd = sprintf('mex %s -lut -lm amd_1.o   amd_2.o   amd_aat.o   amd_control.o   amd_defaults.o   amd_dump.o   amd_global.o   amd_info.o   amd_order.o   amd_post_tree.o   amd_postorder.o   amd_preprocess.o   amd_valid.o     ldl.o   kkt.o   preproc.o   spla.o   cone.o   ecos.o ctrlc.o timer.o   splamm.o   equil.o  ecos_mex.o ecos_bb.o ecos_bb_preproc.o wright_omega.o expcone.o -output "ecos"', d);
        eval(cmd);
    elseif( isunix )
        cmd = sprintf('mex %s -lut -lm -lrt amd_1.o   amd_2.o   amd_aat.o   amd_control.o   amd_defaults.o   amd_dump.o   amd_global.o   amd_info.o   amd_order.o   amd_post_tree.o   amd_postorder.o   amd_preprocess.o   amd_valid.o     ldl.o   kkt.o   preproc.o   spla.o   cone.o   ecos.o ctrlc.o timer.o   splamm.o   equil.o  ecos_mex.o ecos_bb.o ecos_bb_preproc.o wright_omega.o expcone.o -output "ecos"', d);
        eval(cmd);
    end
    fprintf('\t\t\t\t[done]\n');
        
%     fprintf('Copying MEX file...');
%     clear ecos
%     copyfile(['ecos.',mexext], ['../ecos.',mexext], 'f');
%     copyfile( 'ecos.m', '../ecos.m','f');
%     fprintf('\t\t\t[done]\n');
    disp('ecos successfully compiled. Happy solving!');
end

  
%% clean
if( any(strcmpi(what,'clean')) || any(strcmpi(what,'purge')) )
    fprintf('Cleaning up object files...  ');
    if( ispc ), delete('*.obj'); end
    if( isunix), delete('*.o'); end
    fprintf('\t\t[done]\n');
end

%% purge
if( any(strcmpi(what,'purge')) )
    fprintf('Deleting mex file...  ');
    clear ecos
    binfile = ['ecos.',mexext]; 
    if( exist(binfile,'file') )
        delete(['ecos.',mexext]);
    end
    fprintf('\t\t\t[done]\n');
end
