%% take average 
% last = zeros(1,5);
% j=0;
% avg_data = zeros(size(processed_data));
% for i = 1:size(processed_data,1)
%     if (processed_data(i,1)==last(1)&processed_data(i,2)==last(2)&processed_data(i,3)==last(3)&processed_data(i,4)==last(4))
%         last(5) = last(5)+processed_data(i,5);
%         acc_n = acc_n+1;
%     else
%         if (last(1)~=0)
%             j=j+1;
%             avg_data(j,:) = last;
%             avg_data(j,5) = avg_data(j,5)/acc_n;
%         end
%         last = processed_data(i,:);
%         acc_n=1;
%     end
% end
% avg_data = avg_data(1:max(find(avg_data(:,1)>0)),:);

%% get rid of none kept_dose
% kept_dose=[-4,-5,-6,-7,-8];
% n_dose = size(kept_dose,2);
% clean_data = zeros(size(avg_data));
% j=0;
% for i = 1:size(avg_data,1)
%     if (sum(avg_data(i,2)==kept_dose))
%         j=j+1;
%         clean_data(j,:) = avg_data(i,:);
%     end
% end
% clean_data = clean_data(1:max(find(clean_data(:,1)>0)),:);

%% remove ones without n_dose dose
% data = zeros(size(clean_data));
% j=0;
% for i = 1:(size(clean_data,1)-n_dose+1)
%     if (clean_data(i,1)==clean_data(i+n_dose-1,1)&clean_data(i,3)==clean_data(i+n_dose-1,3)&clean_data(i,4)==clean_data(i+n_dose-1,4))
%         data(j+1:j+n_dose,:) = clean_data(i:i+n_dose-1,:);
%         j=j+n_dose;
%     end
% end
% data = data(1:j,:);

%% pick up data from panel pi, cell ci

%     j=0;
%     data_one = zeros(size(data));
%     for i = 1:size(data,1)
%         if (data(i,3)==1&data(i,4)==1)
%             j=j+1;
%             data_one(j,:) = data(i,:);
%         end
%     end
%     data_one = data_one(1:j,:);

%% compute the average reward for each arm, and optimal reward
% mean(data_one(1:n_dose:end,5))
% mean(data_one(2:n_dose:end,5))
% mean(data_one(3:n_dose:end,5))
% 
%     frac=zeros(1,5);
%     mean(data_one(1:5:end,5))
%     mean(data_one(2:5:end,5))
%     mean(data_one(3:5:end,5))
%     mean(data_one(4:5:end,5))
%     mean(data_one(5:5:end,5))
%     best=[];
%     the_sum = zeros(5,1);
%     n_arm = zeros(1,5);y = zeros(1,5);
%     for i = 1:5:size(data_one,1)
%         [~,k] = min(data_one(i:i+4,5));
%         if (k==5)
%             if (binornd(1,0.1)<1)
%                 continue;
%             end
%         end
%         the_sum = the_sum+data_one(i:i+4,5);
%         frac(k)=frac(k)+1;
%         best = [best min(data_one(i:i+4,5))];
%         %random pick an arm to assign
%         arm = randi([1 5],1);
%         n_arm(arm) = n_arm(arm)+1;
%         y(arm) = y(arm)+data_one(i+arm-1,5);
%     end
%     mean(best)
%     the_sum/size(best,2)
%     frac/size(data_one,1)*5
%     y./n_arm
% dif(t) = mean(best)- mean(data_one(5:5:end,5));
% end


%% print SMILES to file
% fileID = fopen('Molecules.csv','w');
% for i = 1:size(SMILES,1)
%     fprintf(fileID,'%i,%s\n',SMILES{i,1},SMILES{i,2});
% end
% fclose(fileID);

%% create features map
% NSC_feature_map_raw = containers.Map('KeyType','double','ValueType','any');
% for i = 2:size(fingerprints,1)
%     if (isnan(fingerprints{i,1}))
%         continue;
%     end
%     NSC_feature_map_raw(fingerprints{i,1}) = fingerprints{i,2};
% end

% % %  %% subsample data
% % % 
% % % %  x=[];
% % % %  y=[];
% % % %  n_arm = zeros(1,n_dose);
% % % % for i = 1:n_dose:size(data_one,1)
% % % %     [~,k] = min(data_one(i:i+n_dose-1,5));
% % % %     arm = randi([1 n_dose],1);
% % % %     n_arm(arm) = n_arm(arm)+1;
% % % %     x(arm, n_arm(arm)) = data_one(i+arm-1,1);
% % % %     y(arm, n_arm(arm)) = data_one(i+arm-1,5);
% % % % end
% % % % [mean(y(1,1:n_arm(1))) mean(y(2,1:n_arm(2))) mean(y(3,1:n_arm(3))) mean(y(4,1:n_arm(4))) mean(y(5,1:n_arm(5)))]

%% dimension reduction

n=size(data_one,1)/n_dose;
target_dimension = 1000;
A=randn(target_dimension,1024);%dimension reduction matrix
big_x = zeros(size(data_one,1)/n_dose,target_dimension);
big_y = zeros(size(data_one,1)/n_dose,n_dose);
k=0;
for i = 1:n_dose:size(data_one,1)
    if (isKey(NSC_feature_map_raw,data_one(i,1)))
            k=k+1;
            s = NSC_feature_map_raw(data_one(i,1));
            int_f = strfind(s,'1');
            feature = sum(A(:,int_f),2);
            big_x(k,:) = feature;
            big_y(k,:) = data_one(i:i+n_dose-1,5)';
    end
end
big_x=big_x(1:k,:);
big_y=big_y(1:k,:);
the_mean = mean(big_x,1);
big_x = big_x-repmat(the_mean,k,1);%centralize x

%[mean(big_y(1:the_split(1))) mean(big_y(the_split(1)+1:the_split(2))) mean(big_y(the_split(2)+1:the_split(3))) mean(big_y(the_split(3)+1:the_split(4))) mean(big_y(the_split(4)+1:the_split(5)))]

% split1 = diff(processed_data(:,1));
% theVar=0;
% total=0;
% best=0;
% nsc_set=[];
% j=1;
% m=zeros(1,5);
% for i =1:size(split1,1)
%     res=[];
%     while(split2(j)<split1(i))
%         res = [res mean(processed_data(split2(j)+1:split2(j+1),5))];
%         j=j+1;
%     end
%     if (size(res,2)==5)
%         nsc_set = [nsc_set processed_data(split1(i),1)];
%         if (res(end)==min(res))
%             if(binornd(1,0.05)<1)
%                 continue;
%             end
%         end
%         if (res(end-1)==min(res))
%             if(binornd(1,0.3)<1)
%                 continue;
%             end
%         end
%         total=total+1;
%         m = m+res;
%         best = best+min(res);
%     end
% end
% size(unique(nsc_set))
% best/total
% m/total

%% remove text
% processed_data = zeros(size(data,1),5);
% nn= size(data,1);
% for i = 2:nn
%     if (mod(i,1000000)==0)
%         i
%     end
%     str_list = split(data{i},',');
%     if (size(str2num(str_list{8})~=0))
%         processed_data(i-1,1) = str2num(str_list{1});
%         processed_data(i-1,2) = str2num(str_list{3});
%         processed_data(i-1,3) = str2num(str_list{6});
%         processed_data(i-1,4) = str2num(str_list{7});
%         processed_data(i-1,5) = str2num(str_list{8});
%     end
% end