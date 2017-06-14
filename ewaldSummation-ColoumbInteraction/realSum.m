% Short-range/Real-space sum contribution

function URealSum = realSum(a1,a2,a3,r,q,N,nReal,L,alpha)

    nL = makePeriodicBox(a1, a2, a3, L, nReal);
    rrCut = (L/2)^2;
    
    URealSum = 0;
    for kk=1:nReal % periodic box #
        for ii=1:N
        for jj=1:N
            dist = norm( r(:,ii) - r(:,jj) + nL(:,kk) );
            if (dist~=0)% && dist^2<=rrCut)
                URealSum = URealSum + q(ii)*q(jj)/dist * erfc(dist*alpha);
            end
        end
        end
    end
    % halve to account for double counting, i.e., ii-jj and jj-ii pair
    % interaction is the same but has been counted twice
    URealSum= URealSum/2;
end