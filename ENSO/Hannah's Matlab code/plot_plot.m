close all
p_k = 21;
p_p = 21;
N_normal = zeros(p_k,p_p);
N_smooth = zeros(p_k,p_p);
% threshold1 = 0.9; % the threshold for absolutly non-chaotic
% threshold2 = 0.1; % the threshold for absolutly chaotic
% for k = 1:p_k
%     for j = 1:p_p
%         if abs(M_normal(k,j))>threshold1
%             N_normal(p_k-k+1,j,:) = zeros(3,1);
%         elseif (threshold1>abs(M_normal(k,j)))&&(abs(M_normal(k,j))>threshold2)
%             N_normal(p_k-k+1,j,:) = [0.5,0.5,0.5];
%         elseif abs(M_normal(k,j))<threshold2
%             N_normal(p_k-k+1,j,:) = ones(3,1);
%         end
%         if abs(M_smooth(k,j))>threshold1
%             N_smooth(p_k-k+1,j,:) = zeros(3,1);
%         elseif (threshold1>abs(M_smooth(k,j)))&&(abs(M_smooth(k,j))>threshold2)
%             N_smooth(p_k-k+1,j,:) = [0.5,0.5,0.5];
%         elseif abs(M_smooth(k,j))<threshold2
%             N_smooth(p_k-k+1,j,:) = ones(3,1);
%         end
%     end
% end
for k = 1:p_k
    for j = 1:p_p
        N_normal(p_k-k+1,j) = M_normal(k,j);
        N_smooth(p_k-k+1,j) = M_smooth(k,j);
    end
end
figure(1)
pcolor(linspace(0,1,21),linspace(1,2,21),M_normal)
colorbar;
xlabel('phi')
ylabel('kappa')
title('Additive and Multiplicative noise with normal time series')
figure(2)
pcolor(linspace(0,1,21),linspace(1,2,21),M_smooth)
colorbar;
xlabel('phi')
ylabel('kappa')
title('Additive and Multiplicative noise with smooth time series')