% Self Interaction term

function USelf = selfInterac(q,eps0,sigma)
    USelf = sum(q.^2) /(4*pi*eps0*sigma*sqrt(2*pi));
end