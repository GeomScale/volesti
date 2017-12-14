function [X,Y] = meshGrid(P, N)
%
%  MESHGRID: Generate X-Y grid for 2D bounded polyhedra. 
%  ======================================================
%  
%  
%  SYNTAX
%  ------
%     
%      [X,Y] = P.meshGrid
%      [X,Y] = P.meshGrid(N)
%      [X,Y] = meshGrid(P,N)
%    
%  
%  DESCRIPTION
%  -----------
%     Generates X-Y grid points for plotting of two-dimensional polyhedra.
%  Supported are only bounded and lower-dimensional polyhedra. The output from this
%  function is consistent with Matlab "meshgrid" function used for plotting
%  functions over 2D polyhedra. The argument to this function is the number of
%  points used to grid the polyhedron, default value is 20.
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron given by V- or                
%            H-representation in dimension 2.         
%            Class: Polyhedron                        
%          N The number of points for gridding. The   
%            value must be positive and               
%            integer-valued.                          
%            Class: double                            
%            Default: 20                              
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          X Matrix with x-coordinates of the grid    
%            values.                                  
%            Class: double                            
%          Y Matrix with y-coordinates of the grid    
%            values.                                  
%            Class: double                            
%              
%  
%  
%  SEE ALSO
%  --------
%     grid,  plot,  fplot
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
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

% local debugging option
PLOTON = 0;

% default grid
if nargin<2
    N = 20;
else
    validate_dimension(N);
end

% use P.forEach() for arrays
error(P.rejectArray());

if ~P.isBounded || P.isEmptySet
    error('Can only grid nonempty and bounded polyhedra.');
end
if P.Dim~=2
    error('The mesh grid can be computed only for 2D polytopes.');
end

% must have vertices
P.minVRep();
V = P.V;

% must have H-rep for faster containement test
P.minHRep();

% Percentage of grid point between each vertex
vs = cell(1,2);
vs{1} = sort(V(:,1)); vs{2} = sort(V(:,2));
% Remove duplicates
for j = 1:2
    vs{j}(abs(diff(vs{j})) < MPTOPTIONS.rel_tol) = [];
end
n = cell(1,2);
for i = 1:2
    n{i} = ceil(N * diff(vs{i}) / (max(V(:,i))-min(V(:,i))));
end

for j = 1:2
    x = vs{j}(1) - (max(V(:,j)) - min(V(:,j))) / N;
    for i = 1:length(vs{j})-1,
        t = linspace(vs{j}(i), vs{j}(i+1), max([2 n{j}(i)]));
        x = [x t(1:end-1)];
    end
    x = [x vs{j}(end)];
    x = [x vs{j}(end) + (max(V(:,j)) - min(V(:,j))) / N];
    xx{j} = x(:);
end

[X,Y] = meshgrid(xx{1},xx{2});

% Compute grid locations of the vertices
gl = zeros(size(V,1),2);
for j = 1:2
    for i = 1:size(V,1)
        xt = abs(xx{j}-V(i,j));
        I = find(xt < MPTOPTIONS.rel_tol);
        if isempty(I)
            [~,I] = min(xt);
        end
        I = I(1);
        
        gl(i,3-j) = I;
    end
end

% the points passed to P.contains() must be column vectors
Z = [X(:), Y(:)]';
II = reshape(P.contains(Z), size(X,1), size(X,2));
% II = zeros(size(X));
% for i = 1:size(X,1)
%   for j = 1:size(X,2)
%     II(i,j) = P.contains([X(i,j);Y(i,j)]);
%   end
% end
VI = zeros(size(X));
for i = 1:size(gl,1)
    II(gl(i,1),gl(i,2)) = 0;
    VI(gl(i,1),gl(i,2)) = 1;
end

if PLOTON
    clf; surf(X,Y,0*X); hold on; plot(P, 'wire', 1); view(2)
end

bnd = zeros(size(X));
for i = 2:size(X,1)-1
    for j = 2:size(X,2)-1
        if VI(i,j) == 1,
            continue;
        end
        if II(i,j) == 1,
            continue; 
        end
        n = 0;
        for ix=-1:1:1
            for iy=-1:1:1
                n = n + II(i+ix,j+iy);
            end
        end
        
        if n > 0,
            bnd(i,j) = 1;
        end
    end
end

% % Order the vertices
% Ikeep = II;
% [a,I] = sort(angle(V*[1;sqrt(-1)]));
% gl = [gl(I,:);gl(I(1),:)];
% bndry = [];
% for i = 1:size(gl,1)-1
%   j = i+1;
%   % Moving from x0 -> x1
%   x0 = [X(gl(i,1),gl(i,2)) Y(gl(i,1),gl(i,2))]';
%   x1 = [X(gl(j,1),gl(j,2)) Y(gl(j,1),gl(j,2))]';
%
%   if PLOTON
%     plot(x0(1),x0(2),'r.','markersize',20);
%     plot(x1(1),x1(2),'b.','markersize',20);
%     hold on;
%   end
%
%   % Deal with the orthogonal cases
%   if gl(i,1) == gl(j,1) || gl(i,2) == gl(j,2)
%     n = max(abs(gl(i,:)-gl(j,:)));
%     t1 = [gl(i,1):gl(j,1)]';
%     t2 = [gl(i,2):gl(j,2)]';
%     l = max([length(t1) length(t2)]);
%     if length(t1) == 1, t1 = t1*ones(l,1); end
%     if length(t2) == 1, t2 = t2*ones(l,1); end
%     t = [t1 t2];
%     bndry = [bndry;t(1:end-1,:)];
%     continue
%   end
%
%   % Compute two directions - one that will move us outside the set and one
%   % inside
%   ang = angle((x1-x0)'*[1;sqrt(-1)]);
%   if ang < 0, ang = ang + 2*pi; end
%
%   if ang < pi/2 % first orthant
%     dOut = [1;0];
%     dIn  = [0;1];
%   elseif ang < pi
%     dOut = [0;1];
%     dIn  = [-1;0];
%   elseif ang < 3*pi/2
%     dOut = [-1;0];
%     dIn  = [0;-1];
%   else
%     dOut = [0;-1];
%     dIn  = [1;0];
%   end
%   dIn = [dIn(2);dIn(1)];
%   dOut = [dOut(2);dOut(1)];
%
%   ii = gl(i,:);
%   while 1
%     bndry(end+1,:) = ii;
%     iIn = ii + dIn';
%     iOut = ii + dOut';
%
%     if norm(iIn-gl(j,:)) == 0
%       break
%     end
%
%     if PLOTON
%       try
%       hIn  = plot(X(iIn(1),iIn(2)),Y(iIn(1),iIn(2)),'c.','markersize',20);
%       hOut = plot(X(iOut(1),iOut(2)),Y(iOut(1),iOut(2)),'k.','markersize',20);
%       catch
%       end
%     end
%
% %     if iIn(1)==0 || iIn(1) > size(X,1) | iIn(2)==0 | iIn(2)>size(X,2)
% %       ii = iOut;
% %     else
%       if II(iIn(1),iIn(2)) == 0
%         ii = iIn;
%       else
%         ii = iOut;
%       end
% %     end
%
%     if PLOTON
%       try
%       hii  = plot(X(ii(1),ii(2)),Y(ii(1),ii(2)),'r.','markersize',20);
%       fprintf('.');
%       delete(hIn); delete(hOut);
%       catch
%       end
%     end
%
%   end
% end

% Move the boundary points onto the boundary
A = P.A;
b = P.b;
R = [eye(2);-eye(2)];
xx = zeros(4,2);
dist = zeros(4,1);
Ikeep = II + bnd + VI;
for ix = 1:size(X,1)
    for iy = 1:size(X,2)
        if bnd(ix,iy) == 0,
            continue;
        end
        
        x = [X(ix,iy);Y(ix,iy)];
        
        % Compute closest point in each direction
        for j = 1:4
            r = R(j,:)';
            pos = find(A*r > MPTOPTIONS.rel_tol);
            t = min((b(pos)-A(pos,:)*x) ./ (A(pos,:)*r));
            if isempty(t)
                t=0;
            end
            xx(j,:) = x'+t*r';
            dist(j) = abs(t);
        end
        
        % Choose closest point
        [~,j] = min(dist);
        xnew = xx(j,:);
        
        % if there are some equalities present, project on the hyperplane
        if size(P.He,1)>0
            nAe = null(P.Ae)';
            xnew = [P.Ae;nAe]\[P.be;nAe*xnew'];
        end
        
        X(ix,iy) = xnew(1);
        Y(ix,iy) = xnew(2);
        
        if PLOTON
            plot([x(1);xnew(1)],[x(2);xnew(2)],'b','linewidth',3);
        end
        
    end
end

X(Ikeep==0) = NaN;
Y(Ikeep==0) = NaN;

% X(:,[1,end])   = [];
% Y(:,[1,end])   = [];
% X([1,end],:)   = [];
% Y([1,end],:)   = [];

end
