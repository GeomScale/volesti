function x = sanitize_inf(x)
% replaces occurences of +/-Inf by +/- MPTOPTIONS.infbound

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

r = MPTOPTIONS.infbound;
x(x==Inf) = r;
x(x==-Inf) = -r;
