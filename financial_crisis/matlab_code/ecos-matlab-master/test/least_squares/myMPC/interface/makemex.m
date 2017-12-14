% FORCES - Fast interior point code generation for multistage problems.
% Copyright (C) 2011-12 Alexander Domahidi [domahidi@control.ee.ethz.ch],
% Automatic Control Laboratory, ETH Zurich.
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

mex -c -O -DUSEMEXPRINTS ../src/myMPC.c
mex -c -O myMPC_mex.c
if( ispc )
    mex myMPC.obj myMPC_mex.obj -output "myMPC"
    delete('*.obj');
elseif( ismac )
    mex myMPC.o myMPC_mex.o -output "myMPC"
    delete('*.o');
else % we're on a linux system
    mex myMPC.o myMPC_mex.o -output "myMPC" -lrt
    delete('*.o');
end
copyfile(['myMPC.',mexext], ['../../myMPC.',mexext], 'f');
copyfile( 'myMPC.m', '../../myMPC.m','f');
