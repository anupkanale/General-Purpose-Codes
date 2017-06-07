% Direct Summation

function UDirect = directSum(a1,a2,a3,r,q,N,nBoxes,L,eps0)

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
    UDirect = UDirect/(4*pi*eps0*2);
end