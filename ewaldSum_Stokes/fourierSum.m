function uFourier = fourierSum(nStokes,rVecM,r_fVec,fVec,xi,L,a1,a2,a3,nImag)
    uFourier = 0;
    
    tau0 = L(1)*L(2)*L(3); % volume of primary cell
    b1 = 2*pi*cross(a2,a3)/L(1);
    b2 = 2*pi*cross(a3,a1)/L(2);
    b3 = 2*pi*cross(a1,a2)/L(3);
    kMat = getP(b1,b2,b3,nImag,L);
    
    for kk=1:(2*nImag+1)^3
        kVec = kMat(:,kk);
        kMag = norm(kVec);
        
        if abs(kMag)>1e-15
            B = 8*pi/kMag^4* (1+kMag^2/(4*xi^2))*(kMag^2*eye(3) - kVec*kVec');
            fFourier = 0; % fourier transformed force vector
            for ii=1:nStokes
                fFourier = fFourier + fVec(:,ii)...
                    * exp(-1j* kVec'*(r_fVec(:,ii)-rVecM) );
            end
            
            coeff = B*exp(-kMag^2/(4*xi^2));
            uFourier = uFourier + real(coeff*fFourier);
        end
    end
    uFourier = uFourier/tau0;
end