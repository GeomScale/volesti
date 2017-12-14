function ok = PWAwelldefined(pwa, Hdom, Kdom, lpsolver, plot)

% 1.) check if the polyhedral partition covers the whole xr-ur space defined
%     by Hdom*[xr;ur] <= Kdom and
% 2.) check if polyhedra overlap 
% 3.) check feasiblity of polyhedra (non-empty interior)

ok = 1;

if ~iscell(Hdom)
    Hdom = {Hdom};
    Kdom = {Kdom};
end
Pdom = polytope;
for ii=1:length(Hdom)
    Pdom = [Pdom polytope(Hdom{ii}, Kdom{ii})];
end

for k=1:length(pwa)
    
    % the k-th pwa partition
    Hi = pwa{k}.Hi;
    Ki = pwa{k}.Ki;
    Pi = Pdom(k);
    
    % dimension of xr-ur space
    dim = size(Hi{1},2);
    
    % find the part of the xr-ur space that is not covered by Hi, Ki
    Options.lpsolver = lpsolver;
    Qdiff = regiondiff(Pi, Pdom);
    
    if isfulldim(Qdiff)
        disp('  Warning: Found gaps in the partition of the x-u space.');
        ok = 0;
        
        % plot the gaps
        if plot > 0,
            if dim == 1
                tg_plotRegions1dim(Qdiff.Hn, Qdiff.Kn);
            elseif dim == 2
                tg_plotRegions(Qdiff.Hn, Qdiff.Kn);
            else
                error('can not plot over this dimension');
            end;
            title('gaps in the partition of the x-u space')
        end;
    end;
    
    % find pairs of overlapping regions
    Qoverl = [];
    for i=1:length(Hi)
        for j=i+1:length(Hi)
            [x, R] = polyinnerball([Hi{i};Hi{j}], [Ki{i};Ki{j}]);
            if R > 1e-6, 
                fprintf('  Warning: regions %i and %i are overlapping.\n', i, j);
                ok = 0;
            end;
        end;
    end;
    
    % check feasibility of polyhedra
    for i=1:length(Ki)
        Pl = polytope(Hi{i}, Ki{i});
        isemptypoly = ~isfulldim(Pl);
        if isemptypoly
            fprintf('  Warning: set %k: regions %i is infeasible.\n', k, i);
        end;
    end;
    
end;
