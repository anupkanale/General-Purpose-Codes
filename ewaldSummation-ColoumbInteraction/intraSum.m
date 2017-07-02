function UIntra = intraSum(r,q,alpha,M,L)

    Nj = 3;
    UIntra = 0.0;
    for jj=1:M
        for kappa=1:Nj
            jjkappa = 3*(jj-1)+kappa;
            
            for lambda=kappa+1:Nj
                jjlambda = 3*(jj-1)+lambda;
                
                dr = r(:,jjlambda) - r(:,jjkappa);
                for index=1:3 % Periodic BC
                    if dr(index)>L/2
                        dr(index) = dr(index) - L;
                    elseif dr(index)<-L/2
                        dr(index) = dr(index) + L;
                    end
                end
                
                dist = norm(dr);
                temp = q(jjkappa)*q(jjlambda)*erf(alpha*dist)/dist;
                UIntra = UIntra + temp;
            end
        end
    end