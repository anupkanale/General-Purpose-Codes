% Short-range/Real-space sum contribution

function URealSum = realSum(a1,a2,a3,r,q,N,nBoxes,L,eps0,sigma)

    nL = makePeriodicBox(a1, a2, a3, L, nBoxes);

    URealSum = 0;
    for kk=1:nBoxes % periodic box #
        for ii=1:N
        for jj=1:N
            dist = norm( r(:,ii) - r(:,jj) + nL(:,kk) );
            if (dist~=0)
                URealSum = URealSum + q(ii)*q(jj)/dist * erfc(dist/(sqrt(2)*sigma));
            end
        end
        end
    end
    URealSum= URealSum/(4*pi*eps0*2);
end