function nL = makePeriodicBox(iVec, jVec, kVec, L, nTerms)
    
    nx = floor(nTerms/2);
    ny = nx; nz = nx;
    nL = zeros(3,nTerms);
    index=0;
    for n1=-nx:nx
        for n2=-ny:ny
            for n3=-nz:nz
                index = index + 1;
                nL(:,index) = n1*L*iVec + n2*L*jVec + n3*L*kVec;
            end
        end
    end
end

% Questions
% 1. nReal and nImag can't take any arbit value, right? because
% can't take individual vectors, have to take entire layer of boxes into
% account.