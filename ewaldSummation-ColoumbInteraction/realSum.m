% Short-range/Real-space sum contribution

function URealSum = realSum(r,q,N,nReal,L,alpha)
    rCut = L/2;
    URealSum = 0;
%     for kk=1:nReal % sum over periodic boxes
    for jj=1:N
    
    
        for ll=jj+1:N
            rTemp = r;
            for index=1:3 % Periodic BC- pulling particles back in box
                if rTemp(index,ll)>rTemp(index,jj)+L/2
                    rTemp(index,ll) = rTemp(index,ll) - L;
                elseif rTemp(index,ll)<rTemp(index,jj)-L/2
                    rTemp(index,ll) = rTemp(index,ll) + L;
                end
            end

%             rjl = rTemp(:,ll) - rTemp(:,jj);
            dist = norm(rTemp(:,ll) - rTemp(:,jj));
            if dist<=rCut % Minimum Image convention
                temp = q(jj)*q(ll)/dist * erfc(dist*alpha);%*heaviside(rCut-dist);
                URealSum = URealSum + temp;
            end
        end
    end
%     end
end


% Notes on Algorithm
% ---------------------------
% 1. real space sum is truncated by rCut, which equals half the primary box
% length
% 2. this ensures that is is computed only in the primary box