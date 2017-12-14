function h = fplot(obj, varargin)
%
%  FPLOT: Plot a single function over a convex set or over an array of convex sets.
%  ================================================================================
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
%     Plot single function over a convex set. If there are more functions attached
%  to a set, then the string name identifies the function to be plotted. If the
%  function is vector valued, i.e. its range is greater than 1, than the first
%  element of the function is plotted by default. For vector valued functions, use
%  the position property to indicate that you want a different element of the
%  function value to plot. Figure properties, such as color, line width, etc, can
%  be specified with "Property" - "Value" pairs.
%  
%  INPUT
%  -----
%     
%        
%          Set    Any object derived from the ConvexSet    
%                 class, e.g. Polyhedron, YSet, ...        
%                 Class: ConvexSet                         
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
%                   ColorMap  Color map given either as a  
%                    M-by-3 matrix, or as a string.        
%                    Default is 'mpt'. Other available     
%                    options are 'hsv', 'hot', 'gray',     
%                    'lightgray', 'bone', 'copper',        
%                    'pink', 'white', 'flag', 'lines',     
%                    'colorcube','vga', 'jet', 'prism',    
%                    'cool', 'autumn', 'spring', 'winter', 
%                    'summer'.                             
%                   ColorOrder  Either 'fixed' for fixed   
%                    ordering of colors, or 'random' for a 
%                    random order. Default is 'fixed'.     
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
 
 
global MPTOPTIONS

%% error checks
error(nargoutchk(0,1,nargout));
if any([obj.Dim]>=3)
    error('Can only plot functions over 2D sets.');
end
if any(~obj.isBounded)
	error('Can only plot bounded polyhedra.');
end
if any(obj.isEmptySet)
	error('Cannot plot function over empty sets.');
end

%% parsing
if mod(numel(varargin), 2)==0
	% P.fplot('opt1', value1, ...)
	start_idx = 1;
	function_name = '';
elseif numel(varargin)>0
	% P.fplot('function_name', 'option', value, ...)
	start_idx = 2;
	function_name = varargin{1};
end
ip = inputParser;
ip.KeepUnmatched = true;
ip.addParamValue('color',  [], @validate_color);
ip.addParamValue('wire',       false,  @(x) islogical(x) || x==1 || x==0);
ip.addParamValue('linestyle',  '-', @validate_linestyle);
ip.addParamValue('linewidth',  1,   @isnumeric);
ip.addParamValue('edgecolor', 'k', @validate_color);
ip.addParamValue('edgealpha', 1, @(x) isnumeric(x) && x>=0 && x<=1);
ip.addParamValue('alpha',      1, @(x) isnumeric(x) && x>=0 && x<=1);
ip.addParamValue('contour',    false, @(x) islogical(x) || x==1 || x==0);
ip.addParamValue('grid',  20,   @isnumeric);
ip.addParamValue('contourGrid',  30,   @isnumeric);
ip.addParamValue('colormap', 'mpt', @(x) (isnumeric(x) && size(x, 2)==3) || ischar(x)); 
ip.addParamValue('colororder', 'fixed', @(x) isequal(x, 'fixed') || isequal(x, 'random'));
ip.addParamValue('showgrid', false, @islogical);
% show_set=true plots the underlying polyhedron
ip.addParamValue('show_set', false, @islogical);
% position=i plots the i-th element f(i) of f=fun(x)
ip.addParamValue('position', 1, @validate_indexset);

% internal parameters
ip.addParamValue('use_hold', true, @islogical);
ip.addParamValue('use_3dview', true, @islogical);
ip.addParamValue('array_index', 1);

ip.parse(varargin{start_idx:end});
options = ip.Results;

%% validation
[function_name, msg] = obj.validateFunctionName(function_name);
error(msg); % the error is only thrown if msg is not empty

fun = obj(1).getFunction(function_name);
if (isa(fun, 'AffFunction') || isa(fun, 'QuadFunction')) && ...
		options.position>fun.R
	error('The position index must be less than %d.', fun.R+1);
end

%% plotting
% hold the plot for the first element of the array
if options.use_hold
	prevHold = ishold;
	if ~ishold,
		newplot;
	end
	hold('on');
end

% plot the array
tic
h = [];
for i=1:numel(obj)
	if toc > MPTOPTIONS.report_period,
		% refresh the plot every 2 seconds
		drawnow;
		tic;
	end
	if numel(obj)>1
		options.array_index = i;
	end
	hi = obj(i).fplot_internal(function_name, options);
	h = [h; hi];
end

if options.use_hold && ~prevHold
	hold('off');
end
if options.use_3dview && any([obj.Dim] >= 2)
	view(3);
	axis tight
end
grid on

if nargout==0
	clear h
end


end
