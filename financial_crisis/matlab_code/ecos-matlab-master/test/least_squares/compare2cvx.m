i=0;
while( 1 )
    data;
    ecos_solver;     ecos_optval
    cvxsocp_solver;  cvx_optval
    
    if( info_.exitflag > 0 )
        assert(strcmp(cvx_status,'Failed'),'ECOS says infeasible, CVX does not');            
    end
    
    if( info_.exitflag == 0 )
        assert(~isempty(strfind(cvx_status,'Solved')),'ECOS says optimal, CVX does not');
        assert( (abs( ecos_optval - cvx_optval ) / abs(cvx_optval)) < 1e-5,'optval differs more than 1e-5'); 
    end
    i=i+1
end