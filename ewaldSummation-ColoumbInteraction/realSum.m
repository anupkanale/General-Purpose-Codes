% Short-range/Real-space sum contribution

function URealSum = realSum(a1,a2,a3,r,q,N,nReal,L,alpha)
    rCut = L/2;
    URealSum = 0;
    nVec = L*getP(a1,a2,a3,nReal);
    
    for kk=1:(2*nReal+1)^3 % sum over periodic boxes
    for jj=1:N
        for ll=jj+1:N
            rjl = r(:,ll) - r(:,jj) + nVec(:,kk);
            for index=1:3 % Periodic BC- pulling particles back in box
                if rjl(index)>rjl(index)+L/2
                    rjl(index) = rjl(index) - L;
                elseif rjl(index)<rjl(index)-L/2
                    rjl(index) = rjl(index) + L;
                end
            end

            dr = norm(rjl);
            if dr~=0 % && distance<=rCut % Minimum Image convention
                temp = q(jj)*q(ll)/dr * erfc(dr*alpha);%*heaviside(rCut-dist);
                URealSum = URealSum + temp;
            end
        end
    end
    end
end