function [isin, inwhich, closest] = isInside(Pn, x0, Options)
%
%  ISINSIDE: Test if a point is contained inside polyhedron in H-representation. 
%  ==============================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      [isin, inwhich, closest] = isInside(Pn, x0, Options)
%      [isin, inwhich, closest] = Pn.isInside(x0, Options)
%    
%  
%  DESCRIPTION
%  -----------
%     Check if the point x0 is contained inside the polyhedron array Pn. The result
%  it the logical statement if x0 in P  and false otherwise. P must be given in
%  H-representation, otherwise an error is thrown. Note that this operation depends
%  on the settings of absolute tolerance that can be changed in Options settings.
%  
%  INPUT
%  -----
%     
%        
%          P                 Polyhedron in H-representation.          
%                            Class: Polyhedron                        
%          x0                A point in the same dimension as         
%                            polyhedron and given as real column      
%                            vector.                                  
%                            Class: double                            
%          Options           A structure with the option settings for 
%                            point location problem.                  
%                            Class: struct                            
%          Options.abs_tol   Absolute tolerance for checking the      
%                            satisfaction of inequalities and         
%                            equalities.                              
%                            Class: double                            
%                            Default: mptopt.abs_tol                  
%          Options.fastbreak If Pn is the polyhedron array, then do a 
%                            quick stop in the consecutive search     
%                            when x0is contained in any of the        
%                            polyhedrons.                             
%                            Class: logical                           
%                            Allowed values:                          
%                                                                     
%                              true                                   
%                              false                                  
%                                                                     
%                            Default: false                           
%                              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          isin    True if x0 in P  and false otherwise.    
%                  Class: logical                           
%                  Allowed values:                          
%                                                           
%                    true                                   
%                    false                                  
%                                                           
%          inwhich If Pn is an array of polyhedra in the    
%                  same dimension, than isin indicates      
%                  which polyhedra from this array contain  
%                  the point x0.                            
%                  Class: double                            
%          closest If Pn is an array of polyhedra in the    
%                  same dimension and none of polyhedra     
%                  contains x0, then the field closest      
%                  indicates which polyhedra has the        
%                  closest distance for x0 to lie in it's   
%                  interior.                                
%                  Class: double                            
%                    
%  
%  
%  SEE ALSO
%  --------
%     contains
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

if nargin<3
	Options = [];
end
if ~isfield(Options, 'fastbreak')
	Options.fastbreak = false;
end
if ~isfield(Options, 'abs_tol')
	Options.abs_tol = MPTOPTIONS.abs_tol;
end

inwhich = [];
closest = [];
for i = 1:length(Pn)
	if any(Pn(i).H_int*[x0; -1] > Options.abs_tol)
		% not in the inequality Hrep
	elseif ~isempty(Pn(i).He_int) && ...
			any(abs(Pn(i).He_int*[x0; -1]) > Options.abs_tol)
		% not in the equality Hrep
	else
		inwhich = [inwhich i];
		if Options.fastbreak
			isin = true;
			return
		end
	end
end

isin = ~isempty(inwhich);
if ~isin && nargout==3
	closest = closestRegion(Pn, x0);
end

end
