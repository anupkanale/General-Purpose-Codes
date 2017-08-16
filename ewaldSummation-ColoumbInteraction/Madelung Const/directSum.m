% direct method

layerList = 5:150;
M = zeros(length(layerList),1);

parfor ll=1:length(layerList)
    layer = layerList(ll);
    M(ll)=0;
    
    for ii=-layer:layer
        for jj=-layer:layer
            for kk=-layer:layer
                dr = (ii^2 + jj^2 + kk^2)^0.5;
                if dr~=0,  M(ll) = M(ll) + (-1)^(ii+jj+kk)/dr; end
            end
        end
    end
    
end

plot(layerList, M, 'o-', 'linewidth', 1.5);