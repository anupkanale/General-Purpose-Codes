% Long-range/Fourier-space sum contribution

function UFourierSum = waveSum(a1,a2,a3,r,q,N,nImag,L,alpha)
    %Reciprocal Lattice vectors
    Vol = L^3;
    b1 = cross(a2,a3)/dot(a1,cross(a2,a3));
    b2 = cross(a3,a1)/dot(a1,cross(a2,a3));
    b3 = cross(a1,a2)/dot(a1,cross(a2,a3));
    pKVec = getk(b1, b2, b3,120,L);
    kMax = 5;
    
    UFourierSum = 0;
    for kk=1:nImag
        kVec = pKVec(:,kk);
        kMag = norm(kVec);
        if (kMag~=0 && kMag^2<kMax^2+2)
            strucFac = 0;
            for ii=1:N
                strucFac = strucFac + q(ii) * exp(2*pi*1j*dot(kVec,r(:,ii)) );
            end
            temp = exp(-(pi*kMag/alpha)^2)* (abs(strucFac)/kMag)^2;
            UFourierSum = UFourierSum + temp;
        end
    end
    UFourierSum = UFourierSum/(2*pi*Vol);
end

function kVector = getk(b1, b2, b3,n,L)
    kVector = [];
    index=1;
    for n1=-n:n
        for n2=-n:n
            for n3=-n:n
                kVector(:,index) = 1/L*(n1*b1+n2*b2+n3*b3)*1e10;
                index=index+1;
            end
        end
    end
end