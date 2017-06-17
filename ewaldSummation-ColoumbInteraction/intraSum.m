function UIntra = intraSum(r,q,alpha,M,L)

    Nj = 3;
    UIntra = 0.0;
    for jj=1:M
        for kappa=1:Nj
            rTemp = r;
            jjkappa = 3*(jj-1)+kappa;
            
            for lambda=kappa+1:Nj
                jjlambda = 3*(jj-1)+lambda;
                
                for index=1:3 % Periodic BC- pulling particles back in box
                    if rTemp(index,jjlambda)>rTemp(index,jjkappa)+L/2
                        rTemp(index,jjlambda) = rTemp(index,jjlambda) - L;
                    elseif rTemp(index,jjlambda)<rTemp(index,jjkappa)-L/2
                        rTemp(index,jjlambda) = rTemp(index,jjlambda) + L;
                    end
                end
                
                dist = norm(rTemp(:,jjlambda) - rTemp(:,jjkappa));
                temp = q(jjkappa)*q(jjlambda)*erf(alpha*dist)/dist;
                UIntra = UIntra + temp;
            end
        end
    end