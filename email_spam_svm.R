# email spam filter 
# dataset from UCI 


#load the two files into R:
dataset=read.csv("data.csv",header=F,sep=";")
names=read.csv("names.csv",header=F,sep=",")

#set data names of the dataset:
names(dataset)=sapply((1:nrow(names)),function(x) as.character(names[x,1]))

#convert the numeric values then we can do a classification problem
dataset$y=as.factor(dataset$y)

#create train, test dataset
set_number=nrow(dataset)
train_prop=0.7
test_prop=0.3

train=sample(set_number,set_number*train_prop)
test=setdiff(1:set_number,train)

train_set=dataset[train,]
test_set=dataset[test,]

#set up packages:
require(e1071)

#create the SVM mode:

svm.model = svm(y~.,data=train_set)
svm.pred = predict(svm.model,test_set[,1:57])

#svm foncusion matrix
svm.matrix=table(pred=svm.pred,true=test_set[,58])
