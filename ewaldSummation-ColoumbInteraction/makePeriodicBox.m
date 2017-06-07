function nL = makePeriodicBox(iVec, jVec, kVec, L, nBoxes)
    
    nx = floor(nBoxes^(1/3)/2);
    ny = nx; nz = nx;
    nL = zeros(3,nBoxes);
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