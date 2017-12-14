function h = pplot(V,varargin)
% PPLOT(V) Plot the matrix V which contains points rowwise
%
% h = pplot(V, fmt, varargin)
%
% If V contains more than three columns, plot the first three
% If V is a vector, then plots one point
%
% set(h,varargin) is applied before return
%

if nargin < 2, 
  fmt = {''}; 
else
  if isstr(varargin)
    fmt = varargin;
  else
    fmt = varargin;
  end
end

if size(V,2) == 1
  V = [V zeros(size(V,1),1)];
end

if size(V,2) == 2
  h = plot(V(:,1),V(:,2), fmt{:});
else
  h = plot3(V(:,1),V(:,2),V(:,3), fmt{:});
end
% if ~isempty(varargin)
%   for i=1:2:length(varargin)
%     set(h, varargin{i}, varargin{i+1});
%   end
% end
