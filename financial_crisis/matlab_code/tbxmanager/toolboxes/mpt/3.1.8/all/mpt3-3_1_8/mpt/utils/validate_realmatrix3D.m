function ts = validate_realmatrix3D(v)
%
% validates real 3D matrix H(:,:,:)
%

if isnumeric(v) && length(size(v))<=3 && all(isfinite(v(:))) && isreal(v)
    ts=true;
else
    error('Input argument must be a real matrix with maximum dimension 3, i.e. H(:,:,:).');
end

end

