% compute OPT by simulation,
function OPT = simulate(mu,H)
    OPT=0;
    trials = 100000;
    [V,D] = eig(H);
    d = size(D,1);
    D_half = real(D.^(1/2));
    for i = 1:trials
        v = V*D_half*randn(d,1)+mu;
        OPT = OPT+max(v);
    end
    OPT = OPT/trials;
end
