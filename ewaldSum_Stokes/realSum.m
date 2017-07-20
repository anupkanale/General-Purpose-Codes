function uReal = realSum(N,rVecM,rVec,fVec,alpha,L,a1,a2,a3,nReal)
    uReal = 0;
    pVec = getP(a1,a2,a3,nReal);

    for nn=1:(2*nReal+1)^3
    for ii=1:N
        dr = rVecM - rVec(:,ii) + pVec(nn);
        drMag = norm(dr);
        if drMag~=0
            A = pi/alpha^(1.5)*(eye(3)*drMag^2 + dr * dr') * phiFunc(pi*drMag^2/alpha, 0.5);
            A = A- 2/sqrt(alpha)*exp(-pi*drMag^2/alpha)*eye(3);
            
            uReal = uReal + A * fVec(:,ii);
        end
    end
    end
end