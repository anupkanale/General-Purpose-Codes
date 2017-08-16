function phiSelf = selfSum(q,N,xi,loc)
    for ii=1:N
        if ii~=loc
            phiSelf = sqrt(xi/pi)*q(ii);
        end
    end
end