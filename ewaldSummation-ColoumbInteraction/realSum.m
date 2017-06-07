% Short-range/Real-space sum contribution

function URealSum = realSum(a1,a2,a3,r,q,N,nBoxes,L,alpha)

    nL = makePeriodicBox(a1, a2, a3, L, nBoxes);

    URealSum = 0;
    for kk=1:nBoxes % periodic box #
        for ii=1:N
        for jj=1:N
            dist = norm( r(:,ii) - r(:,jj) + nL(:,kk) );
            if (dist~=0)
                URealSum = URealSum + q(ii)*q(jj)/dist * erfc(dist*sqrt(alpha));
            end
        end
        end
    end
    % halve to account for double counting, i.e., ii-jj and jj-ii pair
    % interaction is the same but has been counted twice
    URealSum= URealSum/2;
end