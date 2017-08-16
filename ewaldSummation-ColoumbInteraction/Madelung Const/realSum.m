% Short-range/Real-space sum contribution

function phiReal = realSum(a1,a2,a3,r,q,N,nReal,L,xi,loc)
    phiReal = 0;
    nVec = L*getP(a1,a2,a3,nReal);
    
    for kk=1:(2*nReal+1)^3 % sum over periodic boxes
    for ii=1:N
        dr = r(:,loc) - r(:,ii) + nVec(:,kk);
        drMag = norm(dr);
        if drMag>1e-10
            temp = q(ii)*erfc(drMag*sqrt(xi))/drMag;
            phiReal = phiReal + temp;
        end
    end
    end
    phiReal = phiReal*0.5;
end