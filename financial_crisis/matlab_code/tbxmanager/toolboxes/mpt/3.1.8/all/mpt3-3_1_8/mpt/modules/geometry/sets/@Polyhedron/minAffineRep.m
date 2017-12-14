function obj = minAffineRep(obj)
% MINAFFINEREP Compute a minimum representation for the affine set

if obj.hasHRep
    [He, H] = obj.affineHull();
    He = mpt_minAffineRep(He);
    obj.H_int = H;
    obj.He_int = He;
    
elseif obj.hasVRep
    % V-rep
    He = obj.affineHull();
    He = mpt_minAffineRep(He);
    obj.He_int = He;
end

end
