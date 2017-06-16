% Long-range/Fourier-space sum contribution

function UFourierSum = waveSum(a1,a2,a3,r,q,N,nImag,L,alpha)
    %Reciprocal Lattice vectors
    Vol = L^3;
    b1 = cross(a2,a3)/dot(a1,cross(a2,a3));
    b2 = cross(a3,a1)/dot(a1,cross(a2,a3));
    b3 = cross(a1,a2)/dot(a1,cross(a2,a3));
    pKVec = makePeriodicBox(b1, b2, b3, L, nImag);
    kMax = 5;
    
    UFourierSum = 0;
    index = 1;
    for kk=1:nImag
        kVec = pKVec(:,kk);
        kMag = norm(kVec);
        if (kMag~=0 && kMag^2<kMax^2+2)
            strucFac = 0;
            for ii=1:N
                strucFac = strucFac + q(ii) * exp(2*pi*1j*dot(kVec,r(:,ii)) );
            end
            temp(index) = exp(-(pi*kMag/alpha)^2)* (abs(strucFac)/kMag)^2;
            UFourierSum = UFourierSum + temp(index);
            index=index+1;
        end
    end
    UFourierSum = UFourierSum/(2*pi*Vol);
end

% how is the recirocal lattice vector being determined? how are integers
% being obtained as the NISt website says?