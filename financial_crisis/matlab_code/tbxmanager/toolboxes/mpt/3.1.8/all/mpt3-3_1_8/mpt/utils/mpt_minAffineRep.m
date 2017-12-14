function nHe = mpt_minAffineRep(He)
% Compute a minimum representation for the affine set

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

nHe = He;
if size(He,1) == 0
    return
end

r = rank(He, MPTOPTIONS.abs_tol);
if r == size(He, 1)
    % Already minimal
    return
end 

% Choose r linearly independant rows of He
% Note: We don't do any factorizations here, since that can cause fill-in
% of the matrices.
nHe = zeros(r, size(He, 2));
k = 1;
for i = 1:size(He,1)
    if rank([nHe;He(i,:)], MPTOPTIONS.abs_tol) > k-1
        nHe(k,:) = He(i,:);
        k = k+1;
    end
    if k >r, break; end
end

end
