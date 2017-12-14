function smoothLines(h, tf)
%
% Turn on/off line smoothing for all children of h that have the property
%

return

if nargin < 2, tf = true; end

if tf, smooth = 'on'; else smooth = 'off'; end

for i=1:length(h)
  try
    set(h(i),'linesmoothing',smooth);
    hc = get(h(i),'children');
    for j=1:length(hc)
      smoothLines(hc(j));
    end
  catch
  end
end
