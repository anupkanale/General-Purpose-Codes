function uFourier = fourierSum(N,rVecM,r_fVec,fVec,alpha,L,a1,a2,a3,nImag)
    uFourier = 0;
    
    tau0 = L(1)*L(2)*L(3);
    b1 = cross(a2,a3)/tau0;
    b2 = cross(a3,a1)/tau0;
    b3 = cross(a1,a2)/tau0;
    pVec = getP(b1,b2,b3,nImag,ones(3,1));
    
    for kk=1:(2*nImag+1)^3
        kVec = pVec(:,kk);
        kMag = norm(kVec);
        
        if kMag~=0
            expterm = exp(-2*pi*1j*kVec'*rVecM);
            B = pi*alpha^2/tau0*(eye(3)*kMag^2 + kVec*kVec') * phiFunc(pi*alpha*kMag^2, 1);
            fFourier = 0; % fourier transformed force vector
            for ii=1:N
                fFourier = fFourier + fVec(:,ii)...
                    * exp(2*pi*1j*kVec'*r_fVec(:,ii));
            end

            uFourier = uFourier + expterm*B*fFourier;
        end
    end
end