function M = readEspresso(fname)
%===============================================================================
%
% Title:        readEspresso
%                                                             
% Project:      Optimal merging of polyhedra with the same PWA dynamics
%               
% Purpose:      Read Espresso data from output file
%
% Inputs:       fname:  name of Espresso output file
%
% Outputs:      M: markings corresponding to the one output function
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




% read file in
in_file = textread(fname,'%s','delimiter','\n','whitespace','');

% extract number of inputs and outputs
ni = str2num(in_file{1}(4:end));
no = str2num(in_file{2}(4:end));
if no ~= 1, error('file can handle only one output function!'); end;
nmin = str2num(in_file{3}(4:end));

% allocate memory 
M = NaN*(ones(nmin, ni));
for i=1:no, output{i} = []; end;

% read in M
line_offset = 3;
for i = 1:nmin
    % read in terms relating to inputs (minterm)
    for j = 1:ni
        symbol = str2num(in_file{i+line_offset}(j));
        if isempty(symbol), M(i,j) = 0;
        elseif symbol == 0, M(i,j) = -1;
        elseif symbol == 1, M(i,j) = 1; 
        else error('unknown symbol in Espresso file');
        end;
    end  
end;

% prepare data for output
M = M';

return