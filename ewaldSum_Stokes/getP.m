function periodicVec = getP(v1,v2,v3,n,L)
    periodicVec = zeros(3,(2*n+1)^3);
    index=1;
    for n1=-n:n
        for n2=-n:n
            for n3=-n:n
                periodicVec(:,index) = n1*v1*L(1)+n2*v2*L(2) + n3*v3*L(3);
                index=index+1;
            end
        end
    end
end