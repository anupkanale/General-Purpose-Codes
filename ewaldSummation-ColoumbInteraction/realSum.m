% Short-range/Real-space sum contribution

function UReal = realSum(a1,a2,a3,r,q,N,nReal,L,alpha)
    rCut = L/2;
    UReal = 0;
    nVec = L*getn(a1,a2,a3,nReal);
    
    for kk=1:(2*nReal+1)^3 % sum over periodic boxes
    for jj=1:N
        rTemp = r;
        for ll=1:N
            for index=1:3 % Periodic BC- pulling particles back in box
                if rTemp(index,ll)>rTemp(index,jj)+L/2
                    rTemp(index,ll) = rTemp(index,ll) - L;
                elseif rTemp(index,ll)<rTemp(index,jj)-L/2
                    rTemp(index,ll) = rTemp(index,ll) + L;
                end
            end

            rjl = rTemp(:,ll) - rTemp(:,jj) + nVec(:,kk);
            distance = norm(rjl);
            if distance~=0 && distance<=rCut % Minimum Image convention
                temp = q(jj)*q(ll)/distance * erfc(distance*alpha);
                UReal = UReal + temp;
            end
        end
    end
    end
    UReal = UReal/2;
end

function nVector = getn(a1,a2,a3,n)
    nVector = zeros(3,(2*n+1)^3);
    index=1;
    for n1=-n:n
        for n2=-n:n
            for n3=-n:n
                nVector(:,index) = n1*a1+n2*a2+n3*a3;
                index=index+1;
            end
        end
    end
end