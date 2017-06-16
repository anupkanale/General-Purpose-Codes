% Short-range/Real-space sum contribution

function URealSum = realSum(r,q,N,nReal,L,alpha)
    rCut = L/2;
    URealSum = 0;
    for kk=1:nReal % sum over periodic boxes
    for ii=1:N
    for jj=ii+1:N
        for index=1:3 % Periodic BC- pulling particles back in box
            if r(index,jj)>r(index,ii)+L/2
                r(index,jj) = r(index,jj) - L;
            elseif r(index,jj)<r(index,ii)-L/2
                r(index,jj) = r(index,jj) + L;
            end
        end
        rij = r(:,jj) - r(:,ii);
        dist = norm(rij);
        if dist<=rCut % Minimum Image convention
            temp = q(ii)*q(jj)/dist * erfc(dist*alpha);
            URealSum = URealSum + temp;
        end
    end
    end
    end
end

% Notes on Algorithm
% ---------------------------
% 1. real space sum is truncated by rCut, which equals half the primary box
% length
% 2. this ensures that is is computed only in the primary box