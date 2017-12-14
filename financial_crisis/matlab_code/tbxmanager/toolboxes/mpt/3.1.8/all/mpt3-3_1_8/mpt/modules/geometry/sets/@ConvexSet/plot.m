function h = plot(varargin)
%
%  PLOT: Plot the convex set. 
%  ===========================
%  
%  
%  SYNTAX
%  ------
%     
%      h = plot(Set, 'Prop1', value1, 'Prop2', value2)
%      h = Set.plot('Prop1', value1, 'Prop2', value2)
%    
%  
%  DESCRIPTION
%  -----------
%     Plot the general convex set up to dimension three. Figure properties, such as
%  color, line width, etc, can be specified with "Property" - "Value" pairs.
%  
%  INPUT
%  -----
%     
%        
%          Set    Any object derived from the ConvexSet    
%                 class, e.g. Polyhedron, YSet, ...        
%                 Class: ConvexSet                         
%          Prop1  Specification of figure properties.      
%                 Class: char                              
%                 Allowed values:                          
%                                                          
%                   Grid  With how many gridpoints to grid 
%                    the circle/sphere. Default is 40.     
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
 
 
global MPTOPTIONS

% P.plot() -- plots the set (P can be either a single set or an array)
% P.plot('opt1', value1, ...) -- plot with options
% plot(P1, P2, ...) -- plots multiple sets (or arrays thereof)
% plot(P1, 'opt1', value1, P2, 'opt2', value2, ...) -- plots multiple sets
%                                          with different options for each
% h = plot(...) -- returns handle of the plot
%
% All sets must be in dimension <= 3. Empty sets are ignored (generates an
% empty plot). 
%
% Options:
%  'color': color of the function
%     string/double, default=[] (in which case colormap is employed)
%  'wire': if true, the set is plotted in wire frame
%     logical, default=false
%  'linestyle': style of the set's border
%     string, default='-'
%  'linewidth': width of the set's border
%     double, default = 1
%  'edgecolor': color of edges
%     string, default='k'
%  'edgealpha': transparency of edges (1=opaque, 0=complete transparency)
%     double, default=1
%  'alpha': transparency (1=opaque, 0=complete transparency)
%     double, default=1
%  'marker': marker of the plot
%     string, default = 'none'
%  'markerSize': size of the marker
%     double, default=6
%  'colormap': color map to use
%     string or a function handle, default='mpt'
%  'colororder': either 'fixed' or 'random'
%     string, default='fixed'
%  'grid': density of the grid (only used when plotting YSets)
%     numeric, default=40
%  'showindex': if true, displays index of the plotted element (only for
%               bounded polyhedra in 2D)
%     logical, default=false

% Implementation:
% ---------------
%
% ConvexSet/plot is the main (and only) official plotting method for
% ConvexSet and all derived classes. Hence even Polyhedron/plot goes
% via ConvexSet/plot.
%
% ConvexSet/plot is just a dispatcher whose sole purpose is to:
%  1) perform error checks
%  2) parse input arguments
%  3) validate inputs
%  4) call plot_internal() for each element of the array
%
% The actual plotting is done via plot_internal(). Notes:
%  * This internal helper plots just a single set. Does not handle arrays!
%  * No error checks are performed (was already done in ConvexSet/plot)
%  * No imput parsing (done in ConvexSet/plot)
%  * The syntax must be "h = plot_internal(obj, options)":
%  ** must return a vector of handles
%  ** "options" must be a structure generated by inputParser in
%     ConvexSet/plot.  
%  ** all inputs are mandatory
%  ** must always generate the "h" output
%  * If a class defines a custom plotting behavior, then the helper must be
%    declared as a protected method in the class constructor as follows:
% 	    methods (Access=protected)
%           % just include the function prototype here
% 	    	h = plot_internal(obj, options)
% 	    end
%     Then place the actual implementation directly into the class main
%     directory (do not store it into @class/private/ !)
%  * ConvexSet/plot_internal is very general and plots arbitrary convex
%    sets. That's why we have no YSet/plot_internal.
%
% Notes:
% ------
%  * Do NOT introduce plot() methods in subclasses of ConvexSet! It will
%    not work since the method is declared as sealed.
%  * You MUST declare the prototype of plot_internal as a protected
%    method in the class constructor if you want to redefine it (see
%    Polyhedron.m)

%% error checks
error(nargoutchk(0,1,nargout));

%% split input arguments into objects and corresponding options
[objects, options] = parsePlotOptions('ConvexSet', varargin{:});
if numel(objects)==0
	% no objects to plot
	return
end

%% validation
ip = inputParser;
ip.KeepUnmatched = true;
ip.addParamValue('color', [], @validate_color);
ip.addParamValue('wire', false,  @(x) islogical(x) || x==1 || x==0);
ip.addParamValue('linestyle', '-', @validate_linestyle);
ip.addParamValue('linewidth', 1,   @isnumeric);
ip.addParamValue('edgecolor', 'k', @validate_color);
ip.addParamValue('edgealpha', 1, @(x) isnumeric(x) && x>=0 && x<=1);
ip.addParamValue('alpha', 1, @(x) isnumeric(x) && x>=0 && x<=1);
ip.addParamValue('marker', 'none', @validate_marker);
ip.addParamValue('markerSize', 6, @isnumeric);
ip.addParamValue('colormap', 'mpt', @(x) (isnumeric(x) && size(x, 2)==3) || ischar(x)); 
ip.addParamValue('colororder', 'fixed', @(x) isequal(x, 'fixed') || isequal(x, 'random'));
ip.addParamValue('grid', 40, @isnumeric);
ip.addParamValue('showindex', false, @islogical);

% internal: array_index denotes index of the element to plot
ip.addParamValue('array_index', 1, @isnumeric);

dim3plot = false;
for i = 1:numel(objects)
	if any([objects{i}.Dim]>3)
		error('Cannot plot sets in 4D and higher.');
	elseif ~dim3plot && any([objects{i}.Dim]==3)
		dim3plot = true;
	end
end

options_struct = cell(1, numel(options));
for i = 1:numel(options)
	if mod(numel(options{i}), 2)~=0
		error('Options must be in key/value pairs.');
	end
	ip.parse(options{i}{:});
	options_struct{i} = ip.Results;
end

%% plotting

prevHold = ishold;
if ~ishold,
	newplot;
end
hold('on');

% plot each element of each object
tic
h = [];
for i=1:numel(objects)
	if toc > MPTOPTIONS.report_period,
		% refresh the plot every 2 seconds
		drawnow;
		tic;
	end
	for j = 1:numel(objects{i})
		% index of the element in this array
		hij = plot_internal(objects{i}(j), options_struct{i});
		options_struct{i}.array_index = options_struct{i}.array_index + 1;
		h = [h; hij];
	end
end

if ~prevHold
	hold('off');
end
if dim3plot
	view(3);
	axis tight
end
grid on

% do not return handle if not requested
if nargout==0
	clear h
end

axis tight
% enlarge the axis
a = axis;
for i = 1:2:length(a)
	da = a(i+1)-a(i);
	a(i) = a(i)-da*0.05;
	a(i+1) = a(i+1)+da*0.05;
end
axis(a);

end
