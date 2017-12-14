function writeEspresso(fname, M_on, M_off)
%===============================================================================
%
% Title:        writeEspresso
%                                                             
% Project:      Optimal merging of polyhedra with the same PWA dynamics
%               
% Purpose:      Write merging problem as boolean expression into input file
%
% Inputs:       fname:  name of Espresso input file
%               M_on:  markings corresponding to on minterms 
%               M_off: markings corresponding to off minterms
%               both have the format (#hyperplanes x #polyhedra) \in\{-1,-0.5,0.5,1\} 
%
% Remarks:      in the Espresso file, 0-1 means ~AC
%               this is translated into M(i) = -1 0 1
%               This file is intended only for one output function.
%
% Author:       Tobias Geyer 
%                                                                     
% History:      date        ver.    subject                                       
%               2004.02.09  1.0     initial version 
%
% Contact:      Tobias Geyer
%               Automatic Control Laboratory
%               ETH Zentrum, CH-8092 Zurich, Switzerland
%
%               geyer@control.ee.ethz.ch
%
%               Comments and bug reports are highly appreciated
%
%===============================================================================

% (C) 2003-2004 Tobias Geyer, Automatic Control Laboratory, ETH Zurich,
%          geyer@control.ee.ethz.ch

% ---------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------



% rewrite input data as boolean problem
% get Minterms (#minterms x #inputs) \in\{0,1\}
% if find(M_on==0), error('unexpected 0 in M_on'); end;
% if find(M_off==0), error('unexpected 0 in M_off'); end;
M_on = M_on';
M_off = M_off';

ni = size(M_on,2);      % number of inputs
no = 1;                 % number of outputs

% write the preamble of the Espresso file
fid = fopen(fname, 'wt');
fprintf(fid, '.i %i\n', ni);
fprintf(fid, '.o %i\n', no);
fprintf(fid, '.type fr\n');

% if the file is too big, we can use some options:
% -efast:
% 'Stop after the first EXPAND and IRREDUNDANT operations (i.e., do not
% iterate over the solution'
% this is a lot faster and yields almost as many minterms as without
% e.g. 13s instead of 38s, and 12 minterms instead of 11.
%
% -eonset:
% 'Recompute the ON-set before the minimization. Useful when the PLA has a
% large number of product terms (e.g., an exhaustive list of minterms)'
% this is faster and yields almost as many minterms as without
% e.g. 27s instead of 38s, and 12 minterms instead of 11.

% if we want to get the global optimum:
% -Dexact:
% 'Exact minimization algorithm (guarantees minimum number of product terms, 
% and heuristically minimizes number of literals). Potentially expensive.'

% write the on minterms
for i = 1:size(M_on,1)
    for p=1:ni
        if M_on(i,p)==1, fprintf(fid, '1');
        elseif M_on(i,p)==-1, fprintf(fid, '0'); 
        else fprintf(fid, '-'); end;
    end;
    fprintf(fid, ' 1\n');
end;

% write the off minterms
for i = 1:size(M_off,1)
    for p=1:ni
        if M_off(i,p)==1, fprintf(fid, '1');
        elseif M_off(i,p)==-1, fprintf(fid, '0'); 
        else fprintf(fid, '-'); end;
    end;
    fprintf(fid, ' 0\n');
end;

fclose(fid);