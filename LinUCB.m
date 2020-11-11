function [avg_rwd, avg_conf] = LinUCB(X,Y,alpha)
    X = [X ones(size(X,1),1)];
    n_arms = size(Y,2);
    d = size(X,2);
    A=zeros(n_arms,d,d);
    b=zeros(n_arms,d);
    for i =1:n_arms
        A(i,:,:) = eye(d);
    end
    avg_conf = 0;
    avg_rwd = 0;
    theta = zeros(n_arms, d);
    conf = zeros(n_arms,1);
    p_up = zeros(n_arms,1);
    T = size(X,1);
    window_size = 100;
    for i =1:T
        for j = 1:n_arms
            theta(j,:) = squeeze(A(j,:,:))\(b(j,:)');
            conf(j) = alpha*(X(i,:)*(squeeze(A(j,:,:))^(-1))*X(i,:)')^(1/2);
            p_up(j) = theta(j,:)*X(i,:)'+conf(j);
        end
        [~,picked_arm] = max(p_up);
        r=Y(i,picked_arm);
        A(picked_arm,:,:) = squeeze(A(picked_arm,:,:)) + X(i,:)'*X(i,:);
        b(picked_arm,:) = b(picked_arm,:) + r*X(i,:);
        if (i+window_size>T)
            avg_conf = avg_conf + conf(picked_arm);
            avg_rwd = avg_rwd+r;
        end
    end
    avg_conf = avg_conf/window_size;
    avg_rwd = avg_rwd/window_size;
end
