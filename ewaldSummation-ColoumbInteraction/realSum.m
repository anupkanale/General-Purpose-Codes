% Short-range/Real-space sum contribution

% function UReal = realSum(a1,a2,a3,r,q,N,nReal,L,alpha)
%     rCut = L/2; % cut-off distance for min. image conv.
%     nMat = L*getP(a1,a2,a3,nReal); % periodicity vector
% 
%     UReal = 0.0;
%     for kk=1:(2*nReal+1)^3 % sum over periodic boxes
%     nVec = nMat(:,kk);
%     for jj=1:NN
%         for ll=jj+1:N
%             dr = r(:,ll) - r(:,jj) + nVec;
%             for index=1:3 % Periodic BC
%                 if dr(index)>L/2
%                     dr(index) = dr(index) - L;
%                 elseif dr(index)<-L/2
%                     dr(index) = dr(index) + L;
%                 end
%             end
%             
%             drMag = norm(dr);
% %             if drMag~=0
% %             temp = q(jj)/drMag*erfc(drMag*alpha);
%             temp = q(jj)*q(ll)/drMag*erfc(drMag*alpha)*heaviside(rCut-drMag);
%             UReal = UReal + temp;
% %             end
%         end
%     end
%     end
% end

function UReal = realSum(a1,a2,a3,r,q,N,nReal,L,alpha)
    rCut = L/2; % cut-off distance for min. image conv.
    nMat = L*getP(a1,a2,a3,nReal); % periodicity vector

    UReal = 0.0;
    for kk=1:(2*nReal+1)^3 % sum over periodic boxes
    nVec = nMat(:,kk);
    for jj=1:N
            dr = r(:,jj) + nVec;

%             for index=1:3 % Periodic BC
%                 if dr(index)>L/2
%                     dr(index) = dr(index) - L;
%                 elseif dr(index)<-L/2
%                     dr(index) = dr(index) + L;
%                 end
%             end
            
            drMag = norm(dr);
            if drMag~=0
                temp = q(jj)/drMag*erfc(drMag*alpha);
                UReal = UReal + temp;
            end
        end
    end
end