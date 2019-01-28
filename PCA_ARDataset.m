clear all;
close all;

path=pwd;
dataFolder = strcat(pwd,'/AR_database_cropped/test2');
if ~isdir(dataFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', dataFolder);
  uiwait(warndlg(errorMessage));
  return;
end

for l=1:1:13
    filePattern = fullfile(dataFolder, strcat('*-',string(sprintfc('%02d',l)),'.bmp'));
    bmpFiles = dir(filePattern);
    for k = 1:length(bmpFiles)
      baseFileName = bmpFiles(k).name;
      fullFileName = fullfile(dataFolder, baseFileName);
      fprintf(1, 'Now reading %s\n', fullFileName);
      imageArray = imread(fullFileName);
      dataset(:,k+(l-1)*length(bmpFiles))=reshape(imageArray,1,165*120*3);
    end
end


dataset=dataset';
dataset=double(dataset);
classlabel=[ones(100,1);2*ones(100,1);3*ones(100,1);4*ones(100,1);5*ones(100,1);
    6*ones(100,1);7*ones(100,1);8*ones(100,1);9*ones(100,1);10*ones(100,1);11*ones(100,1);
    12*ones(100,1);13*ones(100,1)];
classlabel=double(classlabel);
z=1;
for l=2:1:6
    avg = zeros(13);
    count=zeros(13);
    correct=zeros(13);
    incorrect=zeros(13);
    for k=5:50:1300      
        true_label = classlabel(k);
        data_subset = dataset([1:k-1 k+1:end],:);        
        subclass = classlabel([1:k-1 k+1:end]);

        flag = -1;
        min = Inf;
        for check = 1:13            
            X = data_subset(subclass(:)==check,:);
            Xmean(check,:) = sum(X,1)./length(X(:,1));
            Xmatmean=repmat(Xmean(check,:),length(X(:,1)),1);
            X=X-Xmatmean;
            Xcov = X*(X')/length(X(:,1));
            [V,D] = eigs(Xcov,l);
            V=X'*V;
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
    tavg(z)= sum(avg(:,1))/13;
    z=z+1;
    fprintf('\n Final %f percent correct \n ', tavg(z-1));
    for i = 1:13
        fprintf('Num eigenvectors = %d ; correct %d out of %d, IC %d percent correctly classified = %f %f \n',l,correct(i),count(i), incorrect(i),avg(i)*100),tavg(z-1);
    end
    
    fprintf('\n');
end



