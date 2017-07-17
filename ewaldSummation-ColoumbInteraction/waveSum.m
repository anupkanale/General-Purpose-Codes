% Long-range/Fourier-space sum contribution

function UFourierSum = waveSum(a1,a2,a3,r,q,N,nImag,L,alpha)
    % Reciprocal Lattice vectors
    Vol = L^3;
    b1 = cross(a2,a3)/dot(a1,cross(a2,a3));
    b2 = cross(a3,a1)/dot(a1,cross(a2,a3));
    b3 = cross(a1,a2)/dot(a1,cross(a2,a3));
    kMat = 1/L*getP(b1,b2,b3,nImag);
    kMax = 10;
    
    UFourierSum = 0.0;
    for kk=1:(2*nImag+1)^3
        kVec = kMat(:,kk);
        kMag = norm(kVec);
        if (kMag~=0 && kMag^2<kMax^2 + 2)
            strucFac = 0;
            for ii=1:N
                strucFac = strucFac + q(ii) * exp(2*pi*1j*dot(kVec,r(:,ii)));
            end
            temp = 1/(2*pi*Vol)*exp(-(pi*kMag/alpha)^2)* (abs(strucFac)/kMag)^2;
            UFourierSum = UFourierSum + temp;
        end
    end
end

% function UFourierSum = waveSum(a1,a2,a3,r,q,N,nImag,L,alpha)
%     % Reciprocal Lattice vectors
%     Vol = L^3;
%     b1 = cross(a2,a3)/dot(a1,cross(a2,a3));
%     b2 = cross(a3,a1)/dot(a1,cross(a2,a3));
%     b3 = cross(a1,a2)/dot(a1,cross(a2,a3));
%     kMat = 1/L*getP(b1,b2,b3,nImag);
%     kMax = 10;
%     
%     UFourierSum = 0.0;
%     for kk=1:(2*nImag+1)^3
%         kVec = kMat(:,kk);
%         kMag = norm(kVec);
%         if (kMag~=0 && kMag^2<kMax^2 + 2)
%             strucFac = 0;
%             for ii=1:N
%                 strucFac = strucFac + q(ii) * exp(2*pi*1j*dot(kVec,r(:,ii)));
%             end
%             temp = 1/(pi*Vol)*exp(-(pi*kMag/alpha)^2)* real(strucFac)/kMag^2;
%             UFourierSum = UFourierSum + temp;
%         end
%     end
% end