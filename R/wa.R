`wa` <-
function(...,d.plot=TRUE,env.trans=FALSE,spec.trans=FALSE,diagno=TRUE)

{
par(mfrow=c(1,1))
data<-list(...)
train_set<-data[[1]]
train_env<-data[[2]]
if (length(data)==3)
        test_set<-data[[3]]
                
if (length(data)<=2)
        test_set<-NA
        
if (env.trans=="log10")
    train_env<-log10(train_env)
if (env.trans=="sqrt")
    train_env<-sqrt(train_env)

if (spec.trans=="sqrt")
    {
        train_set<-sqrt(train_set)
        if (length(data)==3)
            {
            test_set<-sqrt(test_set)
            }    
    }

spec_n_train<-NA
n2_train<-NA
spec_n_test<-NA
n2_test<-NA

if (diagno==TRUE)
    {    
        library(vegan)  
        spec_n_train<-specnumber(train_set)
        n2_train<-1/rowSums((train_set/100)^2)
        if (length(data)==3)
            {
                spec_n_test<-specnumber(test_set)
                n2_test<-1/rowSums((test_set/100)^2)
            }
    }
    

col1<-dim(train_set)[2]
row<-dim(train_set)[1]

new<-merge(train_set,test_set,all=TRUE,sort=FALSE)                      # creating test and training_set with the same rows
new<-new[,(1:col1)]
train_set<-new[(1:row),]
test_set<-new[((row+1):dim(new)[1]),]
col<-dim(train_set)[2]


opt<-c(1:col)
opt.cross<-c(1:col)
env_inf.cross<-c(1:row)
env_inf_train<-c(1:row)
env_inf<-c(1:row)

env_infd.cross<-c(1:row)
error_m<-matrix(nrow=1,ncol=8,dimnames=list(("wa-invdes"),c("RMSE","R2","Ave_Bias","Max_Bias","X_R2","X_Ave_Bias","X_Max_Bias","RMSEP")))


opt<-colSums(train_env*train_set)/colSums(train_set)
    env_inf<-colSums(opt*t(train_set),na.rm=TRUE)/rowSums(train_set,na.rm=TRUE)
        fit.lm<-lm(train_env~env_inf)
            opt_d<-fit.lm$coefficients[1]+fit.lm$coefficients[2]*opt
                env_infd_train<-colSums(opt_d*t(train_set),na.rm=TRUE)/rowSums(train_set,na.rm=TRUE)
env_infd<-colSums(opt_d*t(test_set),na.rm=TRUE)/rowSums(test_set,na.rm=TRUE)


RMSE<-sqrt(sum((train_env-env_infd_train)^2)/length(train_env))
R2<-summary(lm(env_infd_train~train_env))$r.squared

    
    
  



############################    cross validation ################################
ord<-order(train_env)
train_env.o<-train_env[ord]
train_set.o<-train_set[ord,]
cross_number<-round(length(train_env)/10,0)


for (i in 1:11)
{
    c1<-(i*cross_number)-(cross_number-1)
    c2<-i*cross_number
    
    if (c1>max(length(train_env)))
         c1<-max(length(train_env))
    if (c2>max(length(train_env)))
         c2<-max(length(train_env))     
    train_set.cross<-train_set.o[-(c1:c2),]
    train_env.cross<-train_env.o[-(c1:c2)]
    test_set.cross<-train_set.o[(c1:c2),]



opt.cross<-colSums(train_env.cross*train_set.cross)/colSums(train_set.cross)
    env_inf.cross<-colSums(opt.cross*t(train_set.cross),na.rm=TRUE)/rowSums(train_set.cross,na.rm=TRUE)
        fit.lm<-lm(train_env.cross~env_inf.cross)
            opt.cross<-fit.lm$coefficients[1]+fit.lm$coefficients[2]*opt.cross

env_infd.cross[c1:c2]<-colSums(opt.cross*t(test_set.cross),na.rm=TRUE)/rowSums(test_set.cross,na.rm=TRUE)

}

######################################## diagnostic plot #######################################################
dif1<-train_env-env_infd_train
dif2<-train_env.o-env_infd.cross


if (length(data)==3)
{

    if (d.plot==TRUE)
         { 
            split.screen(c(1,3))
            split.screen(c(3,1),screen=1)
            split.screen(c(3,1),screen=2)
            screen(4,new = TRUE)  
                plot(train_env,env_infd_train,ylim=c(0,max(train_env)),xlim=c(0,max(train_env)),xlab="observed environmental factor",ylab="inferred env. factor",col="blue")
                abline(0,1,col="red")
            screen(5,new = TRUE)      
                plot(train_env,dif1,xlab="observed environmental factor",ylab="residuals",col="blue")
                abline(0,0,col="red")
            screen(6,new = TRUE)      
                hist(dif1,main="",xlab="error: (inferred - observed)",freq=FALSE,col="darkgreen")
                lines(seq(min(dif1),max(dif1),0.01),dnorm(seq(min(dif1),max(dif1),0.01),mean(dif1),sqrt(var(dif1))),col="red")
            screen(7,new = TRUE)
                plot(train_env.o,env_infd.cross,ylim=c(0,max(train_env.o)),xlim=c(0,max(train_env.o)),main="WA-regression",col.main="blue",font.main=4,cex.main=2,xlab="observed environmental factor",ylab="inferred env. factor - cross val",col="blue")
                abline(0,1,col="red")
            screen(8,new = TRUE)      
                plot(train_env.o,dif2,xlab="observed environmental factor",ylab="residuals - cross val",col="blue")
                abline(0,0,col="red")
            screen(9,new = TRUE)      
                hist(dif2,main="",xlab="cross val error: (inferred - observed)",freq=FALSE,col="darkgreen")
                lines(seq(min(dif2),max(dif2),0.01),dnorm(seq(min(dif2),max(dif2),0.01),mean(dif2),sqrt(var(dif2))),col="red")
                     
            screen(3,new = TRUE)
            x<-c(1:length(env_infd))
             plot(-x~env_infd,axes=FALSE,type="l",col="darkred",ylab=" ",,xlab="reconstructed env.",xlim=c(0,max(train_env)))#
            box()
            axis(1)
            points(env_infd,-x,type="p",pch=19,cex=0.5)
            close.screen (all=TRUE)
         }
    
}    



if (length(data)==2)
{

    if (d.plot==TRUE)
         { 
            split.screen(c(2,1))
            split.screen(c(1,3),screen=1)
            split.screen(c(1,3),screen=2)
            screen(3,new = TRUE)  
                plot(train_env,env_infd_train,ylim=c(0,max(train_env)),xlim=c(0,max(train_env)),xlab="observed environmental factor",ylab="inferred env. factor ",col="blue")
                abline(0,1,col="red")
            screen(4,new = TRUE)      
                plot(train_env,dif1,xlab="observed environmental factor",main="WA-regression",col.main="blue",font.main=4,cex.main=2,ylab="residuals",col="blue")
                abline(0,0,col="red")
            screen(5,new = TRUE)      
                hist(dif1,main="",xlab="error: (inferred - observed)",freq=FALSE,col="darkgreen")
                lines(seq(min(dif1),max(dif1),0.01),dnorm(seq(min(dif1),max(dif1),0.01),mean(dif1),sqrt(var(dif1))),col="red")
            screen(6,new = TRUE)
                plot(train_env.o,env_infd.cross,ylim=c(0,max(train_env.o)),xlim=c(0,max(train_env.o)),xlab="observed environmental factor",ylab="inferred env. factor - cross val",col="blue")
                abline(0,1,col="red")
            screen(7,new = TRUE)      
                plot(train_env.o,dif2,xlab="observed environmental factor",ylab="residuals - cross val",col="blue")
                abline(0,0,col="red")
            screen(8,new = TRUE)      
                hist(dif2,main="",xlab="cross val error: (inferred - observed)",freq=FALSE,col="darkgreen")
                lines(seq(min(dif2),max(dif2),0.01),dnorm(seq(min(dif2),max(dif2),0.01),mean(dif2),sqrt(var(dif2))),col="red")
            
            
    
            close.screen (all=TRUE)        
         }
    
} 


RMSE.X<-sqrt(sum((dif2)^2)/length(train_env.o))
R2.X<-summary(lm(env_infd.cross~train_env.o))$r.squared


error_m[1,1]<-RMSE
error_m[1,2]<-R2
error_m[1,3]<-mean((dif1))
error_m[1,4]<-max(dif1)
error_m[1,5]<-R2.X
error_m[1,6]<-mean((dif2))
error_m[1,7]<-max(dif2)
error_m[1,8]<-RMSE.X
inferred<-"NA"


perfor<-round(error_m,4)

env_infd
inferred<-as.vector(env_infd)


results<-list(spec_n_train,n2_train,opt,env_infd_train,spec_n_test,n2_test,inferred,perfor)
names(results)[[1]]<-"non_z_train"
names(results)[[2]]<-"N2_train"
names(results)[[3]]<-"spec.opt"
names(results)[[4]]<-"inferred.env_train"
names(results)[[5]]<-"non_z_test"
names(results)[[6]]<-"N2_test"
names(results)[[7]]<-"reconstruction"
names(results)[[8]]<-"performance"
results

}

