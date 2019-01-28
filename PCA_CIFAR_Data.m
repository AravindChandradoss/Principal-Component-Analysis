
%4 eigen was best

clear all;
close all;

%Downloads files
if ~exist('cifar-10-batches-mat','dir')
cifar10Dataset = 'cifar-10-matlab';
disp('Downloading 174MB CIFAR-10 dataset...');   
websave([cifar10Dataset,'.tar.gz'],...
    ['https://www.cs.toronto.edu/~kriz/',cifar10Dataset,'.tar.gz']);
gunzip([cifar10Dataset,'.tar.gz'])
delete([cifar10Dataset,'.tar.gz'])
untar([cifar10Dataset,'.tar'])
delete([cifar10Dataset,'.tar'])
end

cifarDataset = 'cifar-10-matlab';
data = load('cifar-10-batches-mat/data_batch_1.mat');

dataset = double(data.data);
classlabel = double(data.labels);
classlabel=classlabel(:,1)+1;
z=1;
for l=5:3:40  
avg = zeros(10);
count=zeros(10);
correct=zeros(10);
incorrect=zeros(10);
fprintf(' For %d eigen vector...running... \n',l);
for k=5:1:60      
    data_subset = dataset([1:k-1 k+1:end],:);        
    subclass = classlabel([1:k-1 k+1:end]);
    true_label = classlabel(k);
    flag = -1;
    min = Inf;
    for check = 1:10            
        X = data_subset(subclass(:)==check,:);
        Xmean(check,:) = sum(X,1)./length(X(:,1));
        Xmatmean=repmat(Xmean(check,:),length(X(:,1)),1);
        X=X-Xmatmean;
        Xcov=X'*X;
        [V,D] = eigs(Xcov,l);
        V=V*1/norm(V,'fro');
        error = norm(dataset(k,:) - (dataset(k,:) - Xmean(check, :))*V*V' - Xmean(check, :),'fro'); 
        if error < min
            min = error;
            flag = check;
        end   

    end
    if flag == true_label
        correct(flag)=correct(flag)+1;
    else
        incorrect(true_label)=incorrect(true_label)+1;
    end
    avg(true_label) = correct(true_label)/(correct(true_label)+incorrect(true_label));
    count(true_label) = count(true_label) +1;
end

tavg(z)= sum(avg(:,1))/10;
z=z+1;
fprintf('\n Final %f percent correct \n ', tavg(z-1));

for i = 1:10
    fprintf('Num eigenvectors = %d ; correct %d out of %d, IC %d percent correctly classified = %f \n',l,correct(i),count(i), incorrect(i),avg(i)*100);
end
fprintf('\n');
end

