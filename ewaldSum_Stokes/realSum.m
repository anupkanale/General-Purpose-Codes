function uReal = realSum(nStokes,rVecM,r_fVec,fVec,alpha,L,a1,a2,a3,nReal)
    uReal = 0;
    pVec = getP(a1,a2,a3,nReal,L);

    for nn=1:(2*nReal+1)^3
    for ii=1:nStokes
        dr = rVecM - r_fVec(:,ii) + pVec(nn);
        
        for index=1:3 % Periodic BC
            if dr(index)>L(index)/2
                dr(index) = dr(index) - L(index);
            elseif dr(index)<-L(index)/2
                dr(index) = dr(index) + L(index);
            end
        end
        drMag = norm(dr);
        
        A = pi/alpha^(1.5)*(drMag^2*eye(3) + dr * dr') * phiFunc(pi*drMag^2/alpha, 0.5);
        A = A- 2/sqrt(alpha)*exp(-pi*drMag^2/alpha)*eye(3);
        uReal = uReal + A * fVec(:,ii);
    end
    end
end