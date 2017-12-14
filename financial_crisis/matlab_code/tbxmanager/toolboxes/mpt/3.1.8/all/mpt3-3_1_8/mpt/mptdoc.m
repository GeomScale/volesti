function mptdoc
%
%  MPTDOC: Display documentation for Multi-Parametric Toolbox in Matlab help
%  =========================================================================
%  browser. 
%  =========
%  
%  
%  SYNTAX
%  ------
%     
%      mptdoc
%    
%  
%  DESCRIPTION
%  -----------
%     The function mptdoc opens the documentation for Multi-Parametric Toolbox
%  inside the Matlab help browser.
%  
%  SEE ALSO
%  --------
%     mpt_init,  mptopt
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Martin Herceg: ETH Zurich
%   mailto:herceg@control.ee.ethz.ch 
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
 
 
v = str2double(strtok(version,'.'));

% locate mptdoc
f=which('mpt.html');

if isempty(f)
    error('The documentation for MPT is not installed on the Matlab path.');
end
if v<8    
    web(f,'-helpbrowser');
else
    % new version
    doc -classic;
end


end
