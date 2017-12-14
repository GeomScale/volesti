function problems = generateData(N)
% Generates testdata.

i = 1;
while( 1 ) 
    data;
    lp;
    cvx_precision best
    if( strcmpi(cvx_status,'solved') || strcmpi(cvx_status,'infeasible') || strcmpi(cvx_status,'unbounded') )
        problems(i).A = A;
        problems(i).b = b;
        problems(i).c = c;
        problems(i).optval = c'*x;
        problems(i).status = cvx_status;
        problems(i).x = x;
        i = i+1;
    end
    if( i > N )
        break;
    end
end
