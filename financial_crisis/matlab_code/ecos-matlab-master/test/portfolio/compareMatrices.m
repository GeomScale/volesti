function ok = compareMatrices(mat1, mat2, threshold)

fprintf('Checking %s and %s...\t', mat1, mat2);

mat1 = dlmread(mat1);
mat2 = dlmread(mat2);

ok = 1;
err = NaN;

% check dimensions
if( size(mat1,1) ~= size(mat2,1) )
    ok = 0;
    status = 'dim1 mismatch';
end
if( size(mat1,2) ~= size(mat2,2) )
    ok = 0;
    status = 'dim2 mismatch';
end

if( ok == 1 )

    % return norm
    err = norm(mat1 - mat2);
    
    if( err <= threshold )
        status = 'OK';
    else
        status = 'failed';
    end
end

if( ~strcmpi(status,'failed') )
    fprintf('%s (%4.2e)\n', status, err);
else
    fprintf(2,'%s (%4.2e)\n', status, err);
end