function [pVec] = makePeriodicBox(iVec, jVec, kVec, L, nBoxes)
    
    nx = floor(nBoxes^(1/3)/2);
    ny = nx; nz = nx;
    pVec = zeros(3,nBoxes);
    index=0;
    for n1=-nx:nx
        for n2=-ny:ny
            for n3=-nz:nz
                index = index + 1;
                pVec(:,index) = n1*L*iVec + n2*L*jVec + n3*L*kVec;
            end
        end
    end
    
%% Questions
%   1. Why is neutral charge condition important?
%   2. Clarify self interaction
%   3. Role of SD of Gaussian cloud