`pom` <-
function(...,d.plot="TRUE")
{
par(mfrow=c(1,1))
data<-list(...)
    train_set<-data[[1]]
    env<-data[[2]]
    if (length(data)==3)
        test_set<-data[[3]]
           
    if (length(data)<=2)
        test_set<-NA
        

#############################################################################################################
#############################################################################################################

pom_int<-function(...,d.plot="TRUE")
{
    data<-list(...)
    train_set<-data[[1]]
    env<-data[[2]]
    if (length(data)==3)
        test_set<-data[[3]]
           
    if (length(data)<=2)
        test_set<-NA

    library(MASS)
    library(mgcv)


    if (length(data)==3)
        {

        col1<-dim(train_set)[2]
        row<-dim(train_set)[1]
        row1<-dim(test_set)[1]
        new<-merge(train_set,test_set,all=TRUE,sort=FALSE)                      # creating test and training_set with the same rows
        new<-new[,(1:col1)]
        train_set<-new[(1:row),]
        test_set<-new[((row+1):dim(new)[1]),]
        col<-dim(train_set)[2]
        train_set<-apply(train_set,2,function(x){ifelse(is.na(x),0,x)})
        test_set<-new[((row+1):dim(new)[1]),]
        test_set<-apply(test_set,2,function(x){ifelse(is.na(x),0,x)})

        }

test_tp<-c(1:3000)
row<-dim(train_set)[2]
col1<-dim(train_set)[1]
if (length(data)==3)
    col2<-dim(test_set)[1]
if (length(data)<=2)
    col2<-2
end_tp<-c(1:col1)
env_mi<-min(env)
env_ma<-max(env)
delta<-(env_ma-env_mi)/3000
for (i in 1:3000) test_tp[i]<-env_mi+i*delta

nspec<-rep(NA,585)
prob_spec<-as.list(nspec)
names(prob_spec)<-names(train_set)
prob_spec.env<-as.list(nspec)
end_tp<-c(1:col1)
tp_all<-matrix(NA,nrow=3000,ncol=col1)

#######################################################

for (i in 1:row)

{    

    fit.glm<-"NA"
    spec.10<-train_set[,i]+10

    spec.pa=ifelse(spec.10>70,7,spec.10)
    spec.pa=ifelse(spec.pa>50,6,spec.pa)
    spec.pa=ifelse(spec.pa>30,5,spec.pa)
    spec.pa=ifelse(spec.pa>20,4,spec.pa)
    spec.pa=ifelse(spec.pa>15,3,spec.pa)
    spec.pa=ifelse(spec.pa>11,2,spec.pa)
    spec.pa=ifelse(spec.pa>10,1,spec.pa)
    spec.pa=ifelse(spec.pa>9,0,spec.pa)
    
    spec.test<-ifelse(spec.pa>0,1,0)

    

if(length(levels(as.factor(spec.pa)))>2)
    
    {
        fit.vglm<-"NA"
        x1<-"NA"
        if (colSums(as.matrix(spec.pa))>=20)
            {
            fit2.polr<-polr(factor(spec.pa)~env+env^2+env^3)
            prob_spec[[i]]<-as.data.frame(fitted(fit2.polr))
            }
    }



}



for (i in 1:584)                                        ################## prediction for environmental parameter #############
    {       
    if(is.null(dim(prob_spec[[i]])[2])==FALSE)
        {
            envfit<-matrix(NA,nrow=3000,ncol=dim(prob_spec[[i]])[2])
            for (k in 1:dim(prob_spec[[i]])[2])
                {
                env1<-env
                fit.gam<-gam(prob_spec[[i]][,k]~s(env1))
                env1<-test_tp
                envfit[,k]<-predict.gam(fit.gam,newdata=as.data.frame(env1),type="response")    
                }
            prob_spec.env[[i]]<-envfit
        }
    }




############################ predicton on the train_set ###############################################


train1_set<-as.matrix(train_set)


for (kn in 1:col1)

{
sample.df<-as.vector(train1_set[kn,])

                            # transformation der Probedaten
    spec.10<-sample.df+10
    spec.10<-as.vector(spec.10)
    spec_pa=as.vector(ifelse(spec.10>70,7,spec.10))
    spec_pa=ifelse(spec_pa>50,6,spec_pa)
    spec_pa=ifelse(spec_pa>30,5,spec_pa)
    spec_pa=ifelse(spec_pa>20,4,spec_pa)
    spec_pa=ifelse(spec_pa>15,3,spec_pa)
    spec_pa=ifelse(spec_pa>11,2,spec_pa)
    spec_pa=ifelse(spec_pa>10,1,spec_pa)
    spec_pa=ifelse(spec_pa>9,0,spec_pa)

sample.pa<-spec_pa


end_tp_try<-rep(1,3000)
for (i in 1:row)
{
    if(is.na(sum((prob_spec.env[[i]])))==FALSE)
        {
           if (dim(prob_spec.env[[i]])[2]>=(sample.pa[i]+1))
            end_tp_try<-end_tp_try*prob_spec.env[[i]][,(sample.pa[i]+1)]
        }
}
tpdat<-cbind(end_tp_try,test_tp)
#plot(end_tp_try~test_tp)
tpdat<-tpdat[order(tpdat[,1]),] 
tp_all[,kn]<-end_tp_try
end_tp[kn]<-tpdat[3000,2]
}


######################################################################### prediction test set #########################################

end_tp_test<-"NA"
tp_all_test<-"NA"
end_tp_test_set<-"NA"

if (length(data)==3)
{
    tp_all_test<-matrix(NA,nrow=3000,ncol=col2+1) 
    train1_set<-as.matrix(test_set)
    end_tp_test_set<-c(1:row1)
    col1<-dim(train1_set)[1]

    for (kn in 1:col1)

    {
        sample.df<-as.vector(train1_set[kn,])

                            # transformation der Probedaten
        spec.10<-sample.df+10
        spec.10<-as.vector(spec.10)
        spec_pa=as.vector(ifelse(spec.10>70,7,spec.10))
        spec_pa=ifelse(spec_pa>50,6,spec_pa)
        spec_pa=ifelse(spec_pa>30,5,spec_pa)
        spec_pa=ifelse(spec_pa>20,4,spec_pa)
        spec_pa=ifelse(spec_pa>15,3,spec_pa)
        spec_pa=ifelse(spec_pa>11,2,spec_pa)
        spec_pa=ifelse(spec_pa>10,1,spec_pa)
        spec_pa=ifelse(spec_pa>9,0,spec_pa)
    
        sample.pa<-spec_pa


        end_tp_test<-rep(1,3000)
        for (i in 1:row)
            {
                if(is.na(sum((prob_spec.env[[i]])))==FALSE)
                    {
                       if (dim(prob_spec.env[[i]])[2]>=(sample.pa[i]+1))
                            end_tp_test<-end_tp_test*prob_spec.env[[i]][,(sample.pa[i]+1)]
                    }
            }
        tpdat_test<-cbind(end_tp_test,test_tp)
        #plot(end_tp_test~test_tp,main=kn)
        tpdat_test<-tpdat_test[order(tpdat_test[,1]),] 
        tp_all_test[,kn]<-end_tp_test
        end_tp_test_set[kn]<-tpdat_test[3000,2]
    }


   

}
results_int<-list(env,test_tp,tp_all,end_tp,tp_all_test,end_tp_test_set)

names(results_int)[[1]]<-"env.test"
names(results_int)[[2]]<-"env.test.pred"
names(results_int)[[3]]<-"env.matrix.test"
names(results_int)[[4]]<-"inf.env.train"
names(results_int)[[5]]<-"env.matrix.train"
names(results_int)[[6]]<-"inf.env.test"
results_int


}


##########################################################################
##########################################################################

par(mfrow=c(1,2))

fit1<-pom_int(train_set,env,test_set)
#if (d.plot==TRUE)
#    {
#        if (length(data)==3)
#            {   row1<-dim(test_set)[1]
                
#                tryp<-fit1[[5]]
#                tryp<-t(as.matrix(tryp))
#                for (k in 1:row1)
#                tryp[k,]<-(tryp[k,]-min(tryp[k,]))/(max(tryp[k,]-min(tryp[k,])))
#                sample<-c(1:row1)
#                image(fit1[[2]],sample,t(tryp)[,row1:1],col=gray.colors(1000,start=1, end=0),axes=FALSE, ylab="sample",xlab="inferred env factor",main="reconstruction")
#                axis(1, at = seq(0.5, 2.5, by = 0.5))
#                box()
#            }
#    }   
#########################################################   cross validation

env_cross<-as.vector(rep(NA,dim(train_set)[1]))
ord<-order(env)
env.od<-env[ord]
train_set.od<-train_set[ord,]
cross_number<-round(length(env.od)/10,0)

    
    for (i in 1:11)
        {
        c1<-(i*cross_number)-(cross_number-1)
        c2<-i*cross_number
        if (c1>max(length(env.od)))
             c1<-max(length(env.od))
        if (c2>max(length(env.od)))
             c2<-max(length(env.od))     
        if (c1<c2)
        {
        train_set.o<-train_set.od[-(c1:c2),]
        env.o<-env.od[-(c1:c2)]
        test_set.o<-as.data.frame(train_set.od[(c1:c2),])
        
        fit.cross<-pom_int(train_set.o,env.o,test_set.o,d.plot="FALSE")
        env_cross[c1:c2]<-as.vector(fit.cross$"inf.env.test")
        }
        }

detach("package:MASS")


#####################################        output                                                               

error_m<-matrix(nrow=1,ncol=8,dimnames=list(("pom"),c("RMSE","R2","Ave_Bias","Max_Bias","X_R2","X_Ave_Bias","X_Max_Bias","RMSEP")))
dif1<-env-fit1[[4]]
dif2<-env-env_cross

error_m[1,1]<-sqrt(sum((dif1)^2)/length(dif1))
error_m[1,2]<-summary(lm(env~fit1[[4]]))$r.squared
error_m[1,3]<-mean(dif1)
error_m[1,4]<-max(dif1)
error_m[1,5]<-summary(lm(env~env_cross))$r.squared
error_m[1,6]<-mean(dif2)
error_m[1,7]<-max(dif2)
error_m[1,8]<-sqrt(sum((dif2)^2)/length(dif2))





###################################    plot                                                      

if (d.plot==TRUE)
{
    
   

    par(mfrow=c(2,3))
           
                
                plot(env,fit1[[4]],ylim=c(0,max(env)),xlim=c(0,max(env)),xlab="observed environmental factor",ylab="inferred env. factor",col="blue")
                abline(0,1,col="red")
           
                plot(env,dif1,xlab="observed environmental factor",ylab="residuals",main="POM regression",col.main="blue",font.main=4,cex.main=2,col="blue")
                abline(0,0,col="red")
          
                hist(dif1,main="",xlab="error: (inferred - observed)",freq=FALSE,col="darkgreen")
               lines(seq(min(dif1),max(dif1),0.01),dnorm(seq(min(dif1),max(dif1),0.01),mean(dif1),sqrt(var(dif1))),col="red")
          
                
                plot(env,env_cross,ylim=c(0,max(env)),xlim=c(0,max(env)),xlab="observed environmental factor",ylab="inferred env. factor - cross val",col="blue")
                abline(0,1,col="red")
           
                plot(env,dif2,xlab="observed environmental factor",ylab="residuals - cross val",col="blue")
                abline(0,0,col="red")
       
                hist(dif2,main="",xlab="cross val error: (inferred - observed)",freq=FALSE,col="darkgreen")
                lines(seq(min(dif2),max(dif2),0.01),dnorm(seq(min(dif2),max(dif2),0.01),mean(dif2),sqrt(var(dif2))),col="red")
         
    		
       
       
      
}

par(mfrow=c(1,1))
results<-list(fit1[[1]],fit1[[2]],fit1[[3]],fit1[[4]],fit1[[5]],fit1[[6]],env_cross,error_m)

names(results)[[1]]<-"env.train"
names(results)[[2]]<-"env.test.pred"
names(results)[[3]]<-"env.matrix.train"
names(results)[[4]]<-"inf.env.train"
names(results)[[5]]<-"env.matrix.test"
names(results)[[6]]<-"reconstruction"
names(results)[[7]]<-"inf.env.cross_train"
names(results)[[8]]<-"performance"
results


}

