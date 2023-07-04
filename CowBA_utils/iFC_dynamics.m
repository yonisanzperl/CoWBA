function [FCD_iFC_triu] = iFC_dynamics(iPC)

    [numAreas,~,numTp] = size(iPC);
    % calculating various FCD measures


    % FCD of upper triangular matrix of iFCs with temporal smoothing (Gus?)
    %%
    % slide adjustment for consistency with FCD => saving as matrix and
    % adding the identity diagonal vlues (taking off +1 on line 11: t2=t+1=numTP-2)
    Isubdiag = find(triu(ones(numAreas),1));
    kk = 1;
    for t = 1:numTp-2
        p1 = squeeze(mean(iPC(:,:,t:t+2),3));
        p1_pattern = p1(Isubdiag);
        for t2 = t+1:numTp-2
            p2 = squeeze(mean(iPC(:,:,t2:t2+2),3));
            p2_pattern = p2(Isubdiag);
            FCD_iFC_triu(kk) = dot(p1_pattern,p2_pattern)/norm(p1_pattern)/norm(p2_pattern);
            %FCD_iFC_triu(t,t2-1) = dot(p1_pattern,p2_pattern)/norm(p1_pattern)/norm(p2_pattern);

            kk = kk+1;
        end
    end
    clear p1 p1_pattern p2 p2_pattern k kk
% %% based on Lopez-Gonzalez et al. 2021
% %FCD (phases)
%     for t=1:numTp
%         patt = iFC(:,:,t);
%         pattern(t,:)=patt(Isubdiag);
%     end
% 
%     npattmax=size(pattern,1);
%     kk3=1;
%     for t=1:npattmax-30
%         p1 = mean(pattern(t:t+30,:));
%         for t2=t+1:npattmax-30
%             p2=mean(pattern(t2:t2+30,:));
%             phfcddata(kk3)=dot(p1,p2)/norm(p1)/norm(p2);
%             matrixcdc(t,t2)=dot(p1,p2)/norm(p1)/norm(p2);
%             matrixcdc(t2,t)=dot(p1,p2)/norm(p1)/norm(p2);
%             kk3=kk3+1;
%         end
%     end 
%     FCD_iFC_triu_slide = matrixcdc;
end