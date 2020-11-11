%clear
% m=10;
% d=1000;
% beta = zeros(m,d);
% for i=1:m
%     beta(i,:) = randn(1,d);
% end
% trueOPT = simulate(zeros(m,1),beta*beta')
% 
n_list = [4000];
n_iter = 20;
estOPT = zeros(1,n_iter);
estPol = zeros(1,n_iter); 
mean_list = zeros(size(n_list));
mean_BM = zeros(size(n_list));
std_list = zeros(size(n_list));
split = zeros(1,m);
for i_n = 1:size(n_list,2)
    n = n_list(i_n);
    for iter = 1:n_iter
        X = randn(m*n,d);
        y=[];
        for i=1:m
            split(i) = n*(i-1)+1;
            y_i = X(split(i):split(i)+n-1,:)*beta(i,:)'+randn(n,1);
            y = [y;y_i];
        end
        estOPT(iter) = est_opt_iso(X,y,split);
        if (n>=d)
            estPol(iter) = est_BM(X,y,split)
        end
    end
    mean_BM(i_n) = mean(estPol);
    mean_list(i_n) = mean(estOPT);
    std_list(i_n) = std(estOPT);
end

figure
set(gcf, 'Position', [100, 1000, 390, 280])
hold on
e=plot([0,n_list(end)],[trueOPT,trueOPT],'-')
e.LineWidth=2;
e=errorbar(n_list,mean_list,std_list);
e.LineWidth = 2;
e=plot(n_list(4:end),mean_BM(4:end));
e.LineWidth = 2;

grid on
legend('ground truth','Iso','Learning');
xlabel('Sample Size')
ylabel('Estimated OPT')
axis([0,n_list(end),46,54])
saveas(gcf,'d=1000_sim','pdf')