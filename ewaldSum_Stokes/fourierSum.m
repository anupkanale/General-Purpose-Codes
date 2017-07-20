function uFourier = fourierSum(N,rVecM,rVec,fVec,alpha,L,a1,a2,a3,nImag)
    uFourier = 0;
    
    tau0 = L^3;
    b1 = cross(a2,a3)/tau0;
    b2 = cross(a3,a1)/tau0;
    b3 = cross(a1,a2)/tau0;
    pVec = getP(b1,b2,b3,nImag);
    
    for kk=1:(2*nImag+1)^3
        kVec = pVec(:,kk);
        kMag = norm(kVec);
        
        B = pi*alpha^2/tau0*(eye(3)*kMag^2 + kVec*kVec) * phiFunc(pi*alpha*kMag^2, 1);
        expterm = exp(-2*pi*1j*dot(kVec,rVecM));
        for ii=1:N
            fFourier = fFourier + fVec(:,ii) * exp(2*pi*1j*dot(kVec,rVec(:,ii)));
        end
        
        uFourier = uFourier + expterm*B*fFourier;
    end
end