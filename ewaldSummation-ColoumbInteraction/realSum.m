% Short-range/Real-space sum contribution

function URealSum = realSum(a1,a2,a3,r,q,N,nReal,L,alpha)
    rCut = L/2;
    URealSum = 0;
    nVec = L*getP(a1,a2,a3,nReal);
    
    for kk=1:(2*nReal+1)^3 % sum over periodic boxes
    for jj=1:N
        for ll=jj+1:N
            dr = r(:,ll) - r(:,jj) + nVec(:,kk);
%             if nReal==0
%                 for index=1:3 % Periodic BC- pulling particles back in box
%                     if dr(index)>L/2
%                         dr(index) = dr(index) - L;
%                     elseif dr(index)<-L/2
%                         dr(index) = dr(index) + L;
%                     end
%                 end
%             end

            drMag = norm(dr);
            if drMag>1e-10
                temp = q(jj)*q(ll)/drMag * erfc(drMag*alpha);
                URealSum = URealSum + temp;
            end
        end
    end
    end
end