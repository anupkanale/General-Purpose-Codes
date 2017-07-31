function hiMat = calcOseenTensor(dr,eta)
    drMag = norm(dr);
    dim = length(dr);
    const = 1/(8*pi*eta);
    
    hiMat = const* (eye(dim)/drMag + dr*dr'/drMag^3);
end