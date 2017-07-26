function hiMat = calcOseenTensor(dr,eta)
    drMag = norm(dr);
    const = 1/(8*pi*eta);
    hiMat = const* (eye(2)/drMag + dr*dr'/drMag^3);
end