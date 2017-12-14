function h = plot(varargin)
%
%  PLOT: Plot the union of convex sets. 
%  =====================================
%  
%  
%  SYNTAX
%  ------
%     
%      h = plot(U, 'Prop1', value1, 'Prop2', value2)
%      h = U.plot('Prop1', value1, 'Prop2', value2)
%      h = plot(U1, 'Prop1', value1, ..., U2, 'Prop2', value2, ...)
%    
%  
%  DESCRIPTION
%  -----------
%     Plot the union of general convex sets up to dimension three. Figure
%  "Value" pairs.
%  
%  INPUT
%  -----
%     
%        
%          U      Union object that contains sets derived  
%                 from the ConvexSet class, e.g.           
%                 Polyhedron, YSet, ...                    
%                 Class: Union                             
%          Prop1  Specification of figure properties.      
%                 Class: char                              
%                 Allowed values:                          
%                                                          
%                   Grid  With how many gridpoints to grid 
%                    the circle/sphere for YSet objects.   
%                    Default is 40.                        
%                   Color  The color of the plot specified 
%                    by real RGB vector or a string name   
%                    of the color (e.g. 'gray');           
%                   Wire  Highlight or not the edges of    
%                    the set. Default is false.            
%                   LineStyle  Specify the type of the     
%                    line to plot edges of the set.        
%                    Accepted values are                   
%                    '-',':','-.','--', and'none'.         
%                   LineWidth  Specify the width of the    
%                    line. Default is 1.                   
%                   Alpha  Transparency of the color. The  
%                    value must be inside [0,1] interval.  
%                    Default is 1.                         
%                   Marker  Type of markings to use.       
%                    Allowed values are ".", "o", "x",     
%                    "+", "*", "s", "d", "v", "", "<",     
%                    ">", "p", "h" or "none". Default is   
%                    "none".                               
%                   MarkerSize  The size of the marker.    
%                    Default is 6.                         
%                   ColorMap  Color map given either as a  
%                    M-by-3 matrix, or as a string.        
%                    Default is 'mpt'. Other available     
%                    options are 'hsv', 'hot', 'gray',     
%                    'lightgray', 'bone', 'copper',        
%                    'pink', 'white', 'flag', 'lines',     
%                    'colorcube', 'vga', 'jet', 'prism',   
%                    'cool', 'autumn', 'spring', 'winter', 
%                    'summer'.                             
%                   ColorOrder  Either 'fixed' for fixed   
%                    ordering of colors, or 'random' for a 
%                    random order. Default is 'fixed'.     
%                   ShowIndex  This option is valid only   
%                    for bounded polyhedra in 2D. If true, 
%                    display an index of the plotted       
%                    element. The default choice is false. 
%                                                          
%                                                          
%          value1 Corresponding value to Prop1.            
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
%     fplot
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
%   mailto:colin.jones@epfl.ch 
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
 
 
error(nargoutchk(0,1,nargout));

h = [];

% split input arguments into objects and corresponding options
[objects, options] = parsePlotOptions('Union', varargin{:});
if numel(objects)==0
	% no objects to plot
	return
end

prevHold = ishold;
if ~ishold, 
    newplot;
end
hold on

% plot each object
h = [];
idx = 1; % index such that each elements is plotted in different color
for i = 1:numel(objects)
	for j = 1:numel(objects{i})
		hj = plot_internal(objects{i}(j), idx, options{i}{:});
		idx = idx + objects{i}(j).Num;
		h = [h; hj];
	end
end

if ~prevHold,
    hold off;
end
if nargout==0
	clear h
end

end
