% Direct Summation

function UDirect = directSum(a1,a2,a3,r,q,N,nBoxes,L)

    nL = makePeriodicBox(a1,a2,a3,L,nBoxes);
    UDirect = 0;
    for kk=1:nBoxes % periodic box #
        for ii=1:N
        for jj=1:N
            dist = norm(r(:,ii)-r(:,jj)+nL(:,kk));
            if( dist~=0)
                UDirect = UDirect + q(ii)*q(jj)/dist;
            end
        end
        end
    end
    % halve to account for double counting, i.e., ii-jj and jj-ii pair
    % interaction is the same but has been counted twice
    UDirect = UDirect/2;
end