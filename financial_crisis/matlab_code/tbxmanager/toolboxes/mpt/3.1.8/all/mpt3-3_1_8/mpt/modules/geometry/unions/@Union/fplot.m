function h = fplot(obj, varargin)
%
%  FPLOT: Plot single function over the sets of the Union object. 
%  ===============================================================
%  
%  
%  SYNTAX
%  ------
%     
%      h = Set.fplot()
%      h = Set.fplot('name', 'Prop1', value1, 'Prop2', value2)
%      h = fplot(Set, 'name', 'Prop1', value1, 'Prop2', value2)
%    
%  
%  DESCRIPTION
%  -----------
%     Plot single function over an union of convex sets. If there are more
%  functions attached to a set, then the string name identifies the function to be
%  plotted. If the function is vector valued, i.e. its range is greater than 1,
%  than the first element of the function is plotted by default. For vector valued
%  functions, use the position property to indicate that you want a different
%  element of the function value to plot. Figure properties, such as color, line
%  width, etc, can be specified with "Property" - "Value" pairs.
%  
%  INPUT
%  -----
%     
%        
%          U      Union object that contains sets derived  
%                 from the ConvexSet class, e.g.           
%                 Polyhedron, YSet, ...                    
%                 Class: Union                             
%          name   If there are more functions attached to  
%                 the set, this string indicates the name  
%                 of the function to plot.                 
%                 Class: char                              
%          Prop1  Specification of figure properties.      
%                 Class: char                              
%                 Allowed values:                          
%                                                          
%                   position  For vector valued functions, 
%                    the position indicates which element  
%                    of the function value to plot.        
%                   Grid  With how many gridpoints to grid 
%                    the circle/sphere. Default is 20.     
%                   Color  The color of the plot specified 
%                    by real RGB vector or a string name   
%                    of the color (e.g. 'gray');           
%                   Wire  Highlight the edges of the set.  
%                    Default is false.                     
%                   LineStyle  Specify the type of the     
%                    line to plot edges of the set.        
%                    Accepted values are                   
%                    '-',':','-.','--', and'none'.         
%                   LineWidth  Specify the width of the    
%                    line. Default is 1.                   
%                   Alpha  Transparency of the color. The  
%                    value must be inside [0,1] interval.  
%                    Default is 1.                         
%                   Contour  Add contour graph. Default is 
%                    false.                                
%                   ContourGrid  With how many grid points 
%                    to plot the contour graph. Default is 
%                    30.                                   
%                   show_set  Plot the domain of the       
%                    function. Default is false.           
%                   showgrid  Show the grid inside the     
%                    set. Default is false.                
%                   colormap  Color map to use given as a  
%                    string or a function handle. Default  
%                    is 'mpt'.                             
%                   colororder  Either 'fixed' or          
%                    'random'. Default is 'fixed'.         
%                                                          
%          value1 Assigns value to Prop1.                  
%                   
%  
%  
%  OUTPUT
%  ------
%     
%        
%          h Handle related to graphics object.       
%            Class: handle                            
%              
%  
%  
%  SEE ALSO
%  --------
%     plot
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
 
 
error(nargoutchk(0,1,nargout));
error(obj.rejectArray());

if iscell(obj.Set)
	% plot each element of the set separately
	
	% hold the plot for the first element of the array
	prevHold = ishold;
	if ~ishold,
		newplot;
	end
	hold('on');

	h = [];
	dim = 0;
	for i = 1:numel(obj.Set)
		% tell ConvexSet/fplot to NOT set hold and axis, since doing so for
		% each element can be slow
		hi = obj.Set{i}.fplot(varargin{:}, 'array_index', i, ...
			'use_hold', false, 'use_3dview', false);
		h = [h; hi];
		dim = max(dim, obj.Set{i}.Dim);
	end
	
	% hold off at the end
	if ~prevHold
		hold('off');
	end
	if dim >= 2
		view(3);
		axis tight
	end

else
	% simpler syntax for arrays
	h = obj.Set.fplot(varargin{:});
end

if nargout==0
	clear h
end

end
