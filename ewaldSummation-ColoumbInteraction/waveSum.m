% Long-range/Fourier-space sum contribution

function UFourierSum = waveSum(a1,a2,a3,r,q,N,nImag,L,alpha)
    %Reciprocal Lattice vectors
    Vol = dot(a1,cross(a2,a3));
    b1 = 2*pi*cross(a2,a3)/Vol;
    b2 = 2*pi*cross(a3,a1)/Vol;
    b3 = 2*pi*cross(a1,a2)/Vol;
    pKVec = makePeriodicBox(b1, b2, b3, L, nImag);
    
    UFourierSum = 0;
    for kk=1:5
        kVec = pKVec(:,kk);
        kMag = norm(kVec);
        if kMag~=0
            strucFac = 0;
            for ii=1:N
                strucFac = strucFac + q(ii) * exp(2*pi*1j*dot(kVec,r(:,ii)) );
            end
            UFourierSum = UFourierSum + exp(-(pi*kMag/alpha)^2)* (abs(strucFac)/kMag)^2;
        end
    end
    UFourierSum = UFourierSum/(2*pi*Vol);
end