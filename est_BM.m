%% estimate the bench mark
function res = estBM(x,y,split)
    m = size(split,2);
    split = [split size(y,1)+1];
    for i = 1:m
        mu(i) = mean(y(split(i):split(i+1)-1));
        y(split(i):split(i+1)-1) = y(split(i):split(i+1)-1)-mu(i);
        beta(:,i) = x(split(i):split(i+1)-1,:)\y(split(i):split(i+1)-1);
        %[Mdl,fit_info] = fitrlinear(x(split(i):split(i+1)-1,:),y(split(i):split(i+1)-1),'Learner','leastsquares','FitBias' ,false,'IterationLimit',1000);
        %beta(:,i) = Mdl.Beta;
        %beta(:,i) = ridge(x(split(i):split(i+1)-1,:), y(split(i):split(i+1)-1), )
    end
    res=0;
    for i = 1:size(x,1)
        res = res + max(x(i,:)*beta+mu);
    end
    res = res/size(x,1);
end