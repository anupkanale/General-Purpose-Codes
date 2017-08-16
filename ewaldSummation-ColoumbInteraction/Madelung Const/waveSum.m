% Long-range/Fourier-space sum contribution

function phiFourier = waveSum(a1,a2,a3,r,q,N,nImag,L,xi,loc)
    %Reciprocal Lattice vectors
    Vol = L^3;
    b1 = 2*pi*cross(a2,a3)/dot(a1,cross(a2,a3));
    b2 = 2*pi*cross(a3,a1)/dot(a1,cross(a2,a3));
    b3 = 2*pi*cross(a1,a2)/dot(a1,cross(a2,a3));
    pMat = 1/L*getP(b1,b2,b3,nImag);
    
    phiFourier = 0;
    for kk=1:(2*nImag+1)^3
        tempi = 0;
        kVec = pMat(:,kk);
        kMag = norm(kVec);
        if kMag~=0
            for ii=1:N
                dr = r(:,loc) - r(:,ii);
                tempi = tempi + q(ii)*exp(-kMag^2/(4*xi) + 1j*kVec'*dr);
            end
        phiFourier = phiFourier + 4*pi/kMag^2*tempi/Vol;
        end
    end
end

% function phiFourier = waveSum(a1,a2,a3,r,q,N,nImag,L,xi,loc)
%     %Reciprocal Lattice vectors
%     Vol = L^3;
%     b1 = 2*pi*cross(a2,a3)/dot(a1,cross(a2,a3));
%     b2 = 2*pi*cross(a3,a1)/dot(a1,cross(a2,a3));
%     b3 = 2*pi*cross(a1,a2)/dot(a1,cross(a2,a3));
%     pMat = 1/L*getP(b1,b2,b3,nImag);
%     
%     phiFourier = 0;
%     for ii=1:N
%         tempi = 0;
%         dr = r(:,loc) - r(:,ii);
%         for kk=1:(2*nImag+1)^3
%             kVec = pMat(:,kk);
%             kMag = norm(kVec);
%             if kMag~=0
%                 tempi = tempi + 4*pi/kMag^2*exp(-kMag^2/(4*xi) + 1j*kVec'*dr);
%             end
%         end
%         
%         phiFourier = phiFourier + q(ii)*tempi/Vol;
%     end
% end