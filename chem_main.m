% n_dose = 5;
% n_iter = 10;
% dim = 1000;
% alpha = 1;
% x = big_x;
% y = -big_y;
% n = 500;
% k = n_dose;
% the_mean = mean(x,1);
% x = x-repmat(the_mean,size(x,1),1);%centralize x
% Sigma = (x'*x)/size(x,1);
% x = x*Sigma^(-1/2);
% X = x;
% Y = y;
% 
% % %% establish ground-truth
% % beta = zeros(dim,n_dose);
% % for i = 1:n_dose
% %     bias(i) = mean(Y(:,i));
% %     beta(:,i) = X\(Y(:,i)-bias(i));
% % end
% % [est_gt_rwd est_gt_arm] = max(X*beta+repmat(bias,size(X,1),1),[],2);
% % est_gt_rwd = mean(est_gt_rwd);
% % gt_rwd = 0;
% % for i =1:size(X,1)
% %     gt_rwd = gt_rwd+Y(i,est_gt_arm(i));
% % end
% % gt_rwd = gt_rwd/size(X,1);
% 
% est_beta = zeros(dim,k);
% estOPT = zeros(1,n_iter);
% learning_rwd = zeros(1,n_iter);
% LinUCB_rwd = zeros(n_iter,2);
% n_list = [100 200 500 1000 1500 2000 4000];
% for i_n = 1:size(n_list,2);
%     n = n_list(i_n);
%     for iter = 1:n_iter
%         X_est = zeros(n*k,dim);
%         Y_est = zeros(n*k,1);
%         split = 1;
%         full_ind = randperm(size(X,1),n*k);
% %         for i =1:k
% %             ind = full_ind((i-1)*n+1:i*n)';
% %             split = [split split(end)+n];
% %             X_est(split(end-1):split(end)-1,:) = X(ind,:);
% %             Y_est(split(end-1):split(end)-1,:) = Y(ind,i);
% %             bias(i) = mean(Y(ind,i));
% %             A = X_est(split(end-1):split(end)-1,:);t
% %             b = (Y(ind,i)-bias(i));
% %             est_beta(:,i) = (A'*A+eye(dim))\(A'*b);
% %         end
% %         split = split(1:end-1);
% %         [~, learning_arm] = max(X*est_beta+repmat(bias,size(X,1),1),[],2);
% %         for i =1:size(X,1)
% %             learning_rwd(iter) = learning_rwd(iter)+Y(i,learning_arm(i));
% %         end
% %         learning_rwd(iter) = learning_rwd(iter)/size(X,1);
% %        if (iter == 1)
%          [LinUCB_rwd((i_n-1)*n_iter+iter, 1) LinUCB_rwd((i_n-1)*n_iter+iter, 2)] = LinUCB(X(full_ind,:),Y(full_ind,:),1);
% %        end
%         %estOPT(iter) = est_opt_iso(X_est,Y_est,split);    
%     end
%     %mean_learning(i_n) = mean(learning_rwd);
%     %mean_list(i_n) = mean(estOPT);
%     %std_list(i_n) = std(estOPT);
%     fprintf('Finished n=%d\n',n);
% end
% for i = 1:7
%     LinUCB_rwd_10_mean(i,:) = mean(LinUCB_rwd_10iter((i-1)*10+1:i*10,:),1); 
% end
figure
hold on
e=plot([0,n_list(end)],[(gt_rwd+est_gt_rwd)/2+100,(gt_rwd+est_gt_rwd)/2+100]);
e.LineWidth=2;
e=errorbar(n_list,LinUCB_rwd_10_mean(:,1)+100,LinUCB_rwd_10_mean(:,2));
e.LineWidth=2;
e=errorbar(n_list,mean_list+100,std_list);
e.LineWidth = 2;

grid on
legend('ground truth','LinUCB', 'Estimated Opt (Alg. 1)');
xlabel('Sample Size')
ylabel('Cancer Growth Inhibition')
axis([0,n_list(end),0,100])
set(gcf, 'Position', [0, 0, 300, 300])
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(gcf,'cancer_d=1000.pdf','-dpdf');
