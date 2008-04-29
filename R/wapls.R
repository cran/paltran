
`wapls` <-function(...,n_comp=4,d.plot=TRUE,plot.comp="RMSEP",env.trans=FALSE,spec.trans=FALSE,diagno=TRUE)
{
par(mfrow=c(1,1))

data<-list(...)
train_set<-data[[1]]
train_env<-data[[2]]
#train_set<-spec.spl$MV[4:588]
#train_env<-env.spl$MV$logtp
#test_set<-viel.df


if (length(data)==3)
        test_set<-data[[3]]
 
if (length(data)<=2)
        test_set<-NA
  
   
col1<-dim(train_set)[2]
row<-dim(train_set)[1]
col<-dim(train_set)[2]
if (length(data)<=2)
    row1<-1
if (length(data)==3)
    row1<-dim(test_set)[1]



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
    




new<-merge(train_set,test_set,all=TRUE,sort=FALSE)                      # creating test and training_set with the same rows
new<-new[,(1:col1)]
train_set<-new[(1:row),]
test_set<-new[((row+1):dim(new)[1]),]
train_set<-apply(train_set,2,function(x){ifelse(is.na(x),0,x)})
test_set<-apply(test_set,2,function(x){ifelse(is.na(x),0,x)})
###########################  declaration of variables #########################
u<-c(1:col1)
numb.c<-c(1:n_comp)
env_inf_train<-matrix(ncol=n_comp,nrow=row)
env_inf<-matrix(0,ncol=n_comp,nrow=row1)
env_infp<-matrix(0,ncol=n_comp,nrow=row1)
RMSE<-c(1:n_comp)
nam<-c("comp1","comp2","comp3","comp4","comp5","comp6","comp7","comp8","comp9","comp10","comp11","comp12","comp13","comp14")
nam1<-nam[1:n_comp]
error<-matrix(0,nrow=n_comp,ncol=10,dimnames=list(nam1,c("RMSE","R2","Ave_Bias","Max_Bias","X_R2","X_Ave_Bias","X_Max_Bias","RMSEP","n_comp","MW")))
r.n1<-matrix(ncol=n_comp,nrow=row)
r.n1.pre<-matrix(ncol=n_comp,nrow=row1)
fac<-c(1:n_comp)
env_infd.cross<-matrix(nrow=row,ncol=n_comp)


train_env.c<-train_env-sum(train_env*rowSums(train_set),na.rm=TRUE)/sum(rowSums(train_set,na.rm=TRUE))
r<-train_env.c


#############

for (run in 1:n_comp)
{

    for(i in 1:col)
        u[i]<-sum(r*train_set[,i],na.rm=TRUE)/sum(train_set[,i],na.rm=TRUE)
    r.n<-as.vector(colSums(u*t(train_set),na.rm=TRUE)/rowSums(train_set,na.rm=TRUE))

    if (run>1)
        {
            v<-sum(((rowSums(train_set)*r.n)*r.n1[,(run-1)]))/sum(rowSums(train_set,na.rm=TRUE))
            r.n<-r.n-(v*r.n1[,run-1])
        }
    
    r.n<-as.vector(r.n)
        z<-sum(r.n*rowSums(train_set))/sum(rowSums(train_set,na.rm=TRUE))
            s2<-sum(((r.n-z)^2)*rowSums(train_set))/sum(rowSums(train_set,na.rm=TRUE))
                r.n1[,run]<-(r.n-z)/sqrt(s2)

    score<-r.n1[,1:run]
    fit.lm<-lm(train_env.c~score,weights=(rowSums(train_set)))

    env_inf_train.c<-fitted(fit.lm)
    env_inf_train[,run]<-fitted(fit.lm)+mean(train_env)
    r<-as.vector((-1)*(train_env.c-env_inf_train.c))

    error[run,1]<-sqrt(sum((train_env-env_inf_train[,run])^2)/length(train_env))
    error[run,3]<-mean((train_env-env_inf_train[,run]))
    error[run,4]<-max(train_env-env_inf_train[,run])
    error[run,2]<-summary(lm(train_env~env_inf_train[,run]))$r.squared
    



############################# prediction for test_set ###############################
if (length(data)==3)
    {
    r.n.pre<-colSums(u*t(test_set),na.rm=TRUE)/rowSums(test_set,na.rm=TRUE)
    if (run>1)
    {
       
        r.n.pre<-r.n.pre-(v*r.n1.pre[,run-1])
    }

    
    r.n1.pre[,run]<-(r.n.pre-z)/sqrt(s2)


    for (k in 1:run)
        env_inf[,run]<-env_inf[,run]+fit.lm$coefficients[k+1]*r.n1.pre[,k]
    
    env_inf[,run]<-env_inf[,run]+fit.lm$coefficients[1]
    env_inf[,run]<-env_inf[,run]+mean(train_env)
    #plot(env_inf[,run])
    }
env_infp<- env_inf
}

#############################     cross validation    ###############################

ord<-order(train_env)
train_env.od<-train_env[ord]
train_set.od<-train_set[ord,]
cross_number<-round(length(train_env)/10,0)

cross_run<-n_comp
    {
    for (i in 1:11)
        {
        c1<-(i*cross_number)-(cross_number-1)
        c2<-i*cross_number
        if (c1>max(length(train_env)))
             c1<-max(length(train_env))
        if (c2>max(length(train_env)))
             c2<-max(length(train_env))     
        if (c1==c2)
            c1<-c1-1
        train_set.o<-train_set.od[-(c1:c2),]
        train_env.o<-train_env.od[-(c1:c2)]
        test_set.o<-train_set.od[(c1:c2),]
        test_set.o<-train_set.od[(c1:c2),]
        row.o<-dim(train_set.o)[1]
        row1.o<-dim(test_set.o)[1]
        col.o<-dim(train_set.o)[2]
        r.n.o<-NA
        u<-rep(NA,col.o)
        r.n1.o<-matrix(0,ncol=cross_run,nrow=row.o)
        r.n1.pre.o<-matrix(0,ncol=cross_run,nrow=row1.o)
        env_inf_train.cross<-matrix(0,ncol=n_comp,nrow=row.o)
        env_inf.o<-matrix(0,ncol=n_comp,nrow=row1.o)
       
        train_env.c<-train_env.o-sum(train_env.o*rowSums(train_set.o),na.rm=TRUE)/sum(rowSums(train_set.o,na.rm=TRUE))

        r<-train_env.c


        

        for (run in 1:cross_run)
        {
        for(i in 1:col.o)
        u[i]<-sum(r*train_set.o[,i],na.rm=TRUE)/sum(train_set.o[,i],na.rm=TRUE)

        r.n.o<-as.vector(colSums(u*t(train_set.o),na.rm=TRUE)/rowSums(train_set.o,na.rm=TRUE))
        if (run>1)
            {
                 v<-sum(((rowSums(train_set.o)*r.n.o)*r.n1.o[,(run-1)]))/sum(rowSums(train_set.o,na.rm=TRUE))
                 r.n.o<-r.n.o-(v*r.n1.o[,run-1])
             }
    
        r.n.o<-as.vector(r.n.o)
        z<-sum(r.n.o*rowSums(train_set.o))/sum(rowSums(train_set.o,na.rm=TRUE))
        s2<-sum(((r.n.o-z)^2)*rowSums(train_set.o))/sum(rowSums(train_set.o,na.rm=TRUE))
        r.n1.o[,run]<-(r.n.o-z)/sqrt(s2)
        score<-r.n1.o[,1:run]
        fit.lm<-lm(train_env.c~score,weights=(rowSums(train_set.o)))
        env_inf_train.c<-fitted(fit.lm)
        env_inf_train.cross[,run]<-fitted(fit.lm)+mean(train_env.o)
        r<-as.vector((-1)*(train_env.c-env_inf_train.c))

        

        r.n.pre.o<-colSums(u*t(test_set.o),na.rm=TRUE)/rowSums(test_set.o,na.rm=TRUE)
        if (run>1)
            {
                r.n.pre.o<-r.n.pre.o-(v*r.n1.pre.o[,run-1])
            }
        
        r.n1.pre.o[,run]<-(r.n.pre.o-z)/sqrt(s2)
              
              
        for (k in 1:run)
         env_inf.o[,run]<-env_inf.o[,run]+fit.lm$coefficients[k+1]*r.n1.pre.o[,k]
        env_infd.cross[c1:c2,run]<-env_inf.o[,run]+fit.lm$coefficients[1]+mean(train_env.o)
        }

        }

        for (i in 1:n_comp)
             {
                inf<-as.vector(env_infd.cross[,i])
                dif2<-train_env.od-inf
                error[i,5]<-summary(lm(inf~train_env.od))$r.squared
                error[i,6]<-mean((dif2))
                error[i,7]<-max((dif2))
                error[i,8]<-sqrt(sum((dif2)^2)/length(train_env.od)) 
            }
        error[,9]<-numb.c
}
#################################     plot   ####################################################
if (d.plot==TRUE)
    {
        if (plot.comp=="RMSEP")
            {
            error1<-error[order(error[,8]),]
            comp<-error1[1,9]
            inf<-(env_inf_train[,comp])
            inf.cross<-as.vector(env_infd.cross[,comp])
            dif1<-train_env-inf
            dif2<-train_env-inf.cross

            if (length(data)==2)
                {    
                    split.screen(c(2,1))
                    split.screen(c(1,3),screen=1)
                    split.screen(c(1,3),screen=2)
                   screen(3,new = TRUE)  
                
                    plot(train_env,inf,ylim=c(0,max(train_env)),xlim=c(0,max(train_env)),xlab="observed environmental factor",ylab=paste("inferred env. factor c",comp),col="blue")
                        abline(0,1,col="red")
                
                    screen(4,new = TRUE)      
                    plot(train_env,dif1,xlab="observed environmental factor",ylab=paste("residuals c",comp),main="WA-PLS regression",col.main="blue",font.main=4,cex.main=2,col="blue")
                        abline(0,0,col="red")
                
                    screen(5,new = TRUE)      
                    hist(dif1,main="",xlab="error: (inferred - observed)",freq=FALSE,col="darkgreen")
                        lines(seq(min(dif1),max(dif1),0.01),dnorm(seq(min(dif1),max(dif1),0.01),mean(dif1),sqrt(var(dif1))),col="red")
                
                    screen(6,new = TRUE)
                    plot(train_env,inf.cross,ylim=c(0,max(train_env)),xlim=c(0,max(train_env)),xlab="observed environmental factor",ylab=paste("inferred env. factor c",comp," - cross val"),col="blue")
                        abline(0,1,col="red")
                    
                    screen(7,new = TRUE)      
                    plot(train_env,dif2,xlab="observed environmental factor",ylab=paste("residuals c",comp," - cross val"),col="blue")
                        abline(0,0,col="red")
                    
                    screen(8,new = TRUE)      
                    hist(dif2,main="",xlab="cross val error: (inferred - observed)",freq=FALSE,col="darkgreen")
                        lines(seq(min(dif2),max(dif2),0.01),dnorm(seq(min(dif2),max(dif2),0.01),mean(dif2),sqrt(var(dif2))),col="red")
         
                    close.screen (all=TRUE)    
                } 
            if (length(data)==3)
                {
                    split.screen(c(1,3))
                    split.screen(c(3,1),screen=1)
                    split.screen(c(3,1),screen=2)
                    
                    screen(4,new = TRUE)  
                    plot(train_env,inf,ylim=c(0,max(train_env)),xlim=c(0,max(train_env)),xlab="observed environmental factor",ylab=paste("inferred env. factor c",comp),col="blue")
                        abline(0,1,col="red")
                    
                    screen(5,new = TRUE)      
                    plot(train_env,dif1,xlab="observed environmental factor",ylab=paste("residuals c",comp),col="blue")
                        abline(0,0,col="red")
                    
                    screen(6,new = TRUE)      
                    hist(dif1,main="",xlab="error: (inferred - observed)",freq=FALSE,col="darkgreen")
                        lines(seq(min(dif1),max(dif1),0.01),dnorm(seq(min(dif1),max(dif1),0.01),mean(dif1),sqrt(var(dif1))),col="red")
                    
                    screen(7,new = TRUE)
                    plot(train_env,inf.cross,ylim=c(0,max(train_env)),xlim=c(0,max(train_env)),main="WA-PLS regression",col.main="blue",font.main=4,cex.main=2,xlab="observed environmental factor",ylab=paste("inferred env. factor c",comp," - cross val"),col="blue")
                        abline(0,1,col="red")
                    
                    screen(8,new = TRUE)      
                    plot(train_env,dif2,xlab="observed environmental factor",ylab=paste("residuals c",comp," - cross val"),col="blue")
                        abline(0,0,col="red")
                    
                    screen(9,new = TRUE)      
                    hist(dif2,main="",xlab="cross val error: (inferred - observed)",freq=FALSE,col="darkgreen")
                        lines(seq(min(dif2),max(dif2),0.01),dnorm(seq(min(dif2),max(dif2),0.01),mean(dif2),sqrt(var(dif2))),col="red")

                    screen(3,new = TRUE)
                    inf_test<-env_inf[,comp]
                    x<-c(1:length(inf_test))
                    plot(-x~inf_test,axes=FALSE,type="l",col="darkred",ylab=" ",xlab=paste("reconstructed env. wapls comp ",comp),xlim=c(0,max(train_env)))#
                    box()
                    axis(1)
                    points(inf_test,-x,type="p",pch=19,cex=0.5)
                    close.screen (all=TRUE)
                }
            }
}

############################## output ##############################

results<-list(spec_n_train,n2_train,score,env_inf_train,spec_n_test,n2_test,env_infp,error)
names(results)[[1]]<-"non_z_train"
names(results)[[2]]<-"N2_train"
names(results)[[3]]<-"scores"
names(results)[[4]]<-"inferred.env_train"
names(results)[[5]]<-"non_z_test"
names(results)[[6]]<-"N2_test"
names(results)[[7]]<-"reconstruction"
names(results)[[8]]<-"performance"
results

}