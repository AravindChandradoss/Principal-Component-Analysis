clear all

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
    for l = 1:length(bmpFiles)
      baseFileName = bmpFiles(l).name;
      fullFileName = fullfile(dataFolder, baseFileName);
      fprintf(1, 'Now reading %s\n', fullFileName);
      imageArray = imread(fullFileName);
      dataset(:,l+(l-1)*length(bmpFiles))=reshape(imageArray,1,165*120*3);
    end
end


dataset=dataset';
dataset=double(dataset);
classlabel=[ones(100,1);2*ones(100,1);3*ones(100,1);4*ones(100,1);5*ones(100,1);
    6*ones(100,1);7*ones(100,1);8*ones(100,1);9*ones(100,1);10*ones(100,1);11*ones(100,1);
    12*ones(100,1);13*ones(100,1)];
classlabel=double(classlabel);
z=1;
avg = zeros(13);
count=zeros(13);
correct=zeros(13);
incorrect=zeros(13);


for l=5:2:15
    data_subset = dataset([1:l-1 l+1:end],:);
    sub_class = classlabel([1:l-1 l+1:end]);
    true_label = classlabel(l);
    flag = -1;
    max = - Inf; 
    for check = 1:1            
        X = data_subset(sub_class(:)==check,:);
        Xmean1(check,:) = sum(X,1)./length(X(:,1));
        Xpic= Xmean1(check,:)-(sum(Xmean1(check,:))./length(Xmean1(check,:)));
        Ypic= dataset(l,:)-(sum(dataset(l,:))./length(dataset(1,:)));
        
        corr=(Xpic*Ypic')/(std(Xpic)*std(Ypic));
        
        if corr > max
            max = corr;
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
for i = 1:13
    fprintf(' correct %d out of %d, IC %d percent correctly classified = %f \n',correct(i),count(i), incorrect(i),avg(i)*100);
end
fprintf('\n');


