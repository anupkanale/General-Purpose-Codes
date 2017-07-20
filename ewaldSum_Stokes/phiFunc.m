function phi = phiFunc(x, option)
    if option==0.5
        phi = exp(-x)/x + erfc(sqrt(x))/(2*x);
    elseif option==1
        phi = exp(-x)*(1+x)/x^2;
    else
        disp('Invalid option')
        phi = 0;
    end
end