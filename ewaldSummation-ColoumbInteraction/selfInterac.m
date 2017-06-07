% Self Interaction term

function USelf = selfInterac(q,alpha)
    USelf = sum(q.^2) * sqrt(alpha/pi);
end