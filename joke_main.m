dim = 2000;
k=10;
% filename = 'C:\Users\Jun\Documents\dataset\jester-data-1.xls';
% M1 = xlsread(filename);
% filename = 'C:\Users\Jun\Documents\dataset\jester-data-2.xls';
% M2 = xlsread(filename);
% raw_data = [M1;M2];
% 
% %% take top k jobs as arms, remove people who did not rate at least one of the top k jokes.
% M = raw_data;
% rate_per_joke = sum(M<99,1);
% rate_per_joke = rate_per_joke(2:end);
% [~,I] = sort(rate_per_joke,'descend');
% topK = I(1:k);
% has_missing = find(sum(M(:,topK+1)>10,2)>0);
% M(has_missing,:)=[];
% % 
% %% replace the missing value by the user's mean
% %compute mean
% the_mean = (sum(M(:,2:end),2)-(100-M(:,1))*99)./M(:,1);
% for i = 1:size(M,1)
%     ind = find(M(i,:)>10);
%     M(i,ind) = the_mean(i);
% end
% M = M(:,2:end);
% % 
% %% create Y, remove Y from M
% joke_feature = 90;
% Y = M(:,topK);
% M = M(:,I(k+1:k+joke_feature));
% 
% 
% %% create matrix X,y 
% % X = zeros(size(M,1),joke_feature^2);
% % for i = 1:size(M,1)
% %     x = M(i,:)'*M(i,:);
% %     x = x(:);
% %     X(i,:) = x';
% % end
% % 
% 
% %% create matrix X,y sigmoid version
% X = M*randn(joke_feature,dim)/sqrt(dim);
% X = 1./(1+exp(-X));
% 
% %% random projection
% % X = X*randn(joke_feature^2,dim)/joke_feature;
% 
% %% centralize data
% X = X-repmat(mean(X,1),size(M,1),1);
% 
% %% scale the data
% S = (X'*X)/size(X,1);
% X = X*S^(-1/2);
% fprintf('Preprocessing finished');
% %% establish ground-truth
% beta = zeros(dim,k);
% for i = 1:k
%     bias(i) = mean(Y(:,i));
%     beta(:,i) = X\(Y(:,i)-bias(i));
% end
% [est_gt_rwd est_gt_arm] = max(X*beta+repmat(bias,size(X,1),1),[],2);
% est_gt_rwd = mean(est_gt_rwd);
% gt_rwd = 0;
% for i =1:size(X,1)
%     gt_rwd = gt_rwd+Y(i,est_gt_arm(i));
% end
% gt_rwd = gt_rwd/size(M,1);
% save(['10joke_d=' num2str(dim) '_sig.mat'],'X','Y','beta','bias','gt_rwd');
% % 
% fprintf('Ground truth established');

%% scale the data again
% S = (X'*X)/size(X,1);
% X = X*S^(-1/2);
%% split into k parts
n_list = [100 200 500 1000 1500 2000 2500 3000];
n_iter = 10;
mean_list = zeros(size(n_list));
std_list = zeros(size(n_list));
mean_BM = zeros(size(n_list));
est_beta = zeros(dim,k);
for i_n = 1:size(n_list,2);
    n = n_list(i_n);
    estOPT = zeros(1,n_iter);
    learning_rwd = zeros(1,n_iter);
    for iter = 1:n_iter
%         X_est = zeros(n*k,dim);
%         Y_est = zeros(n*k,1);
%         split = 1;
        full_ind = randperm(size(X,1),n*k);
%         for i =1:k
%             ind = full_ind((i-1)*n+1:i*n)';
%             split = [split split(end)+n];
%             X_est(split(end-1):split(end)-1,:) = X(ind,:);
%             Y_est(split(end-1):split(end)-1,:) = Y(ind,i);
%             bias(i) = mean(Y(ind,i));
%             est_beta(:,i) = X_est(split(end-1):split(end)-1,:)\(Y(ind,i)-bias(i));
%         end
%         split = split(1:end-1);
%         [~, learning_arm] = max(X*est_beta+repmat(bias,size(X,1),1),[],2);
%         for i =1:size(X,1)
%             learning_rwd(iter) = learning_rwd(iter)+Y(i,learning_arm(i));
%         end
%         learning_rwd(iter) = learning_rwd(iter)/size(X,1);

        %estOPT(iter) = est_opt_iso(X_est,Y_est,split);
        LinUCB(X())
        
    end
    mean_BM(i_n) = mean(learning_rwd);
%    mean_list(i_n) = mean(estOPT);
%    std_list(i_n) = std(estOPT);
    fprintf('Finished n=%d\n',n);
end
%save(['d=' num2str(dim) '_sig.mat'],'mean_list','std_list','mean_BM','gt_rwd');
% random-projection with sigmoid
% sample with replacement

% scaled_mean_list = mean_list/4+2.5;
% scaled_std_list = std_list/4;
% scaled_gt_rwd = gt_rwd/4+2.5;
%% create plot
% 
% figure
% set(gcf, 'Position', [100, 1000, 390, 280])
% hold on
% e=plot([0,n_list(end)],[scaled_gt_rwd,scaled_gt_rwd]);
% e.LineWidth=2;
% e=plot(n_list,scaled_mean_BM,'s-','MarkerSize',2)
% e.LineWidth=2;
% e=errorbar(n_list,scaled_mean_list,scaled_std_list);
% e.LineWidth = 2;


% load('d=100_sig.mat')
% e=plot([0,n_list(end)],[gt_rwd,gt_rwd],'-')
% e.LineWidth=1;
% e=errorbar(n_list,mean_list,std_list);
% e.LineWidth = 1;
% load('d=200_sig.mat')
% e=plot([0,n_list(end)],[gt_rwd,gt_rwd],'-')
% e.LineWidth=1;
% e=errorbar(n_list,mean_list,std_list);
% e.LineWidth = 1;

grid on
%legend('ground truth-50','Iso-50','ground truth-100','Iso-100','ground truth-200','Iso-200');
legend('ground truth','Perf. of Learned Policy','Estimated Opt (Alg. 1)');
%title('Estimating OPT')
xlabel('Sample Size')
ylabel('Reward')
axis([0,n_list(end),0,5])
saveas(gcf,'d=2000_sig','pdf')
