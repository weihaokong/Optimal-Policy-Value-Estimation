%% Estimate the expected reward of the optimal policy
%  X is a matrix whose rows correspond to each context 
%  y is a vector contains rewards
%  split is a vector contains the start index of each arm's record

% est_opt_reward_1 is for the isotropic covariance setting
function res = est_opt_iso(X,y,split)
    m = size(split,2);
    split = [split size(y,1)+1];
    H = zeros(m);
    G = triu(X*X',1);
    for i = 1:m
        mu(i) = mean(y(split(i):split(i+1)-1));
        y(split(i):split(i+1)-1) = y(split(i):split(i+1)-1)-mu(i);
    end
    mu=mu';
    for i = 1:m
        for j = i:m
            H(i,j) = y(split(i):split(i+1)-1)'*G(split(i):split(i+1)-1,split(j):split(j+1)-1)*y(split(j):split(j+1)-1);
            if (i==j)
                bot = nchoosek(split(i+1)- split(i)-1,2);
            else
                bot = ((split(i+1)- split(i)-1)*(split(j+1)- split(j)));
            end
            H(i,j) = H(i,j)/bot;
            H(j,i) = H(i,j);
        end
    end
    res = simulate(mu,H);
end
