wa <-
function(...,d.plot=TRUE,env.trans=FALSE,spec.trans=FALSE,diagno=TRUE,val=c("none","n.cross","loo","boot"),run=10,scale=FALSE,seed=1,out=TRUE,desh.meth=c("class","inverse"),drop.non.sig=FALSE,min.occ=1,nfold=10)
{
    
    data<-list(...)
    data_l<-length(data)
    train_set<-as.matrix(data[[1]]) 
    train_set<-train_set[,colSums(train_set)!=0]
    train_env<-as.matrix(data[[2]])
    if (data_l==3)
        test_set<-data[[3]]
    if (data_l<=2)
        test_set<-NA
    
    if (missing(val)) 
        {
        val<- "none"
        val1<-"none"
        }
    if (missing(desh.meth)) 
        desh.meth<-"inverse"
    
    
    if (min.occ>1)
        {
            min.occ.taxa<-apply(train_set,2,function(x) if(sum(ifelse(x>0,1,0))>=min.occ) x<-1 else x<-0)
            train_set<-train_set[,min.occ.taxa==1]
        }

    
    
    
    
    if (drop.non.sig=="TRUE")
        {
            library(mgcv)
            col<-dim(train_set)[2]
            sig.ts<-NA
            for (i in 1:col)
                {
                    env<-train_env[,1]
                    if (sum(ifelse(train_set[,i]>0,1,0))>5)
                        {
                            fit.gam<-gam(train_set[,i]~s(env))    
                            k<-summary(fit.gam)[[8]]
                            if (k<0.05)
                                sig.ts[i]<-1 else sig.ts[i]<-0
                        } else    sig.ts[i]<-0 
                }
            train_set<-train_set[,sig.ts==1]
            train_set1<-train_set[rowSums(train_set)!=0,]
            train_env<-as.matrix(train_env[rowSums(train_set)!=0,])
            train_set<-train_set1
        }
  
   
   
    if (scale==TRUE)
        {
            if (length(data)==3)
            {   
                test_set<-data[[3]]
                if (dim(test_set)[1]>1)
                    test_set<-as.data.frame(t(apply(test_set,1,function(x)x*100/sum(x,na.rm=TRUE)))) 
            }
            train_set<-as.data.frame(t(apply(train_set,1,function(x)x*100/sum(x,na.rm=TRUE))))
        }            

    
    
    if (env.trans=="log10")
        train_env<-log10(train_env)
    if (env.trans=="sqrt")
        train_env<-sqrt(train_env)

    if (spec.trans=="sqrt")
        {
            train_set<-sqrt(train_set)
            if (length(data)==3)
                test_set<-sqrt(test_set)        
        }
    spec_n_train<-NA
    n2_train<-NA
    spec_n_test<-NA
    n2_test<-NA
    RMSEP<-NA
    n_test.train<-NA
    
    if (diagno==TRUE)
        {    
            library(vegan)  
            spec_n_train<-specnumber(train_set)
            n2_train<-1/rowSums((train_set/100)^2)
            if (length(data)==3)
                { 
                   if (dim(test_set)[1]>1)
                        {
                        spec_n_test<-specnumber(test_set) 
                        n2_test<-1/rowSums((test_set/100)^2)
                        k<-colnames(test_set) %in% colnames(train_set)
                        k<-colnames(test_set)[k]
                        test_set<-test_set[,k]                 
                        n_test.train<-specnumber(test_set)
                        }
                }
        }
    
###########################################   wa              ####################################################

wa1<-function(...)
{
    data<-list(...)
    train_set<-data[[1]]
    train_env<-data[[2]]
    dat<-length(data)
    opt<-NA
    if (dat==3)
     {        
        test_set<-data[[3]]
        test_set<-as.matrix(test_set)
     }
    sum1<-rowSums(train_set,na.rm=TRUE)
    t.train_set<-t(train_set)
    opt<-colSums(apply(train_set,2,function(x) train_env*x))/colSums(train_set)           
    env_inf_train<-colSums(apply(t.train_set,2,function(x) opt*x),na.rm=TRUE)/sum1
       
    if (desh.meth=="class")
        {
            fit.lm<-lm(env_inf_train~train_env)
            opt_d<-(opt-fit.lm$coefficients[1])/fit.lm$coefficients[2]
            env_infd_train<-colSums(opt_d*t.train_set,na.rm=TRUE)/sum1
        }
    if (desh.meth=="inverse")
        {
            fit.lm<-lm(train_env~env_inf_train)
            opt_d<-fit.lm$coefficients[1]+fit.lm$coefficients[2]*opt
            env_infd_train<-colSums(opt_d*t.train_set,na.rm=TRUE)/sum1
        }
    inferred<-NA           
    if (dat==3)   
        {
        if (dim(test_set)[1]>=2)
            {
            k<-colnames(test_set) %in% names(opt_d)
            k<-colnames(test_set)[k]
            test_set<-test_set[,k]                                                       
            opt_d<-opt_d[k]
            inferred<-colSums(opt_d*t(test_set),na.rm=TRUE)/rowSums(test_set)
            }
        if (dim(test_set)[1]<2)
            {
            k<-colnames(test_set) %in% names(opt_d)
            k<-colnames(test_set)[k]
            test_set<-test_set[,k]                                                       
            opt_d<-opt_d[k]
            inferred<-sum(opt_d*t(test_set),na.rm=TRUE)/sum(test_set)
            }
        
        }
    dif<-env_infd_train-train_env
    RMSE<-sqrt(mean(dif^2))       
    R2<-summary(fit.lm)$r.squared
    mean_er<-mean((dif))
    max_er<-max(abs(dif))
    results<-list(opt,mean_er,max_er,inferred,R2,RMSE,env_infd_train,opt_d)
    results    
}  # end of wa1

    
###########################        LOO            ##############################  

loo.wa<-function(train_set,train_env)
    {
    dim_ts<-nrow(train_set)
    loo1<-c(1:dim_ts)
    for (i in seq_len(dim_ts))
        {
        train_set.c<-train_set[-i,]
        train_set.c<-train_set.c[,colSums(train_set.c)!=0]
        test_set.c<-train_set[i,]
        k<-names(test_set.c)[names(test_set.c) %in% colnames(train_set.c)]
        test_set.c<-test_set.c[names(test_set.c) %in% colnames(train_set.c)]
        train_set.c<-train_set.c[,k]
        train_env.c<-train_env[-i]
        opt<-wa1(train_set.c,train_env.c)[[8]]
        
        loo1[i]<-sum(opt*test_set.c)/sum(test_set.c)
        }
    names(loo1)<-row.names(train_set)
    dif<-loo1-train_env
    max_er.c<-max(abs(dif))
    mean_er.c<-mean(dif)
    RMSEP<-sqrt(mean(dif^2))
    R2.X<-summary(lm(loo1~train_env))$r.squared
    error<-list(R2.X,mean_er.c,max_er.c,RMSEP,loo1)
    error
    }
 
loo_pred.wa<-function(train_set,train_env,test_set)
    {
    dim_ts<-nrow(train_set)
    loo1<-c(1:dim_ts)
    loo_pred<-matrix(ncol=dim_ts,nrow=dim(test_set)[1])
    sum2<-rowSums(test_set) 
    test_set<-as.matrix(test_set)
    for (i in seq_len(dim_ts))
        {
        train_set.c<-train_set[-i,]
        train_set.c<-train_set.c[,colSums(train_set.c)!=0]
        test_set.c<-train_set[i,]
        k<-names(test_set.c)[names(test_set.c) %in% colnames(train_set.c)]
        test_set.c<-test_set.c[names(test_set.c) %in% colnames(train_set.c)]
        train_set.c<-train_set.c[,k]
        train_env.c<-train_env[-i]
        #opt<-wa1(train_set.c,train_env.c)[[8]]        
        #loo1[i]<-sum(opt*test_set.c)/sum(test_set.c)
        #loo_pred[,i]<-wa1(train_set.c,train_env.c,test_set)[[4]]
        mod<-wa1(train_set.c,train_env.c) 
        opt<-mod[[8]]
        loo1[i]<-sum(opt*test_set.c)/sum(test_set.c)
        k<-colnames(test_set) %in% names(opt)
        k<-colnames(test_set)[k]
        test_set.opt<-test_set[,k]                                                       #Skalierung auf 100%???
        opt_d<-opt[k]
        if(dim(test_set)[1]>1)
            loo_pred[,i]<-colSums(opt_d*t(test_set.opt),na.rm=TRUE)/rowSums(test_set.opt)
        if(dim(test_set)[1]==1)
            loo_pred[,i]<-sum(opt_d*t(test_set.opt),na.rm=TRUE)/sum(test_set.opt)
      }
    names(loo1)<-row.names(train_set)
    dif<-loo1-train_env
    max_er.c<-max(abs(dif))
    mean_er.c<-mean(dif)
    RMSEP<-sqrt(mean(dif^2))
    loo_inf<-apply(loo_pred,1,mean)
    loo_sd<-apply(loo_pred,1,sd)
    names(loo_inf)<-row.names(test_set)
    names(loo_sd)<-row.names(test_set)
    row.names(loo_pred)<-row.names(test_set)
    colnames(loo_pred)<-c(1:dim_ts)
    R2.X<-summary(lm(loo1~train_env))$r.squared
    error<-list(R2.X,mean_er.c,max_er.c,RMSEP,loo1,loo_inf,loo_sd,loo_pred)
    error
    } 
 
 
############################    cross validation ################################

tencross.wa<-function(train_set,train_env,run1=run,n.fold=nfold)
    {
    max_er.c<-NA
    mean_er.c<-NA
    RMSEP<-NA
    R2.c<-NA
    dim_ts<-dim(train_set)[1]
    loo1<-matrix(nrow=dim_ts,ncol=run1)
    loo1.b<-matrix(nrow=dim_ts,ncol=run1)
    loo1.k<-matrix(nrow=dim_ts,ncol=run1)
    c.cross<-NA
    cross_number<-round((dim_ts/n.fold),0)
    set.seed(seed) 
    for (i in seq_len(n.fold))
        c.cross[i]<-(i-1)*cross_number
    c.cross[n.fold+1]<-max(length(train_env))
                
    for (r in seq_len(run1))
    {                  
        k<-sample((1:dim(train_set)[1]))
        train_set.b<-train_set[k,]
        train_env.b<-train_env[k,]
        loo1.b[,r]<-train_env.b
        loo1.k[,r]<-k
        for (i in seq_len(n.fold))
        {
            c1<-c.cross[i]+1
            c2<-c.cross[i+1]
            train_set.c<-train_set.b[-(c1:c2),]
            train_env.c<-train_env.b[-(c1:c2)]
            test_set.c<-train_set.b[(c1:c2),]        
            train_set.c<-train_set.c[,colSums(train_set.c)!=0]
            loo1[c1:c2,r]<-wa1(train_set.c,train_env.c,test_set.c)[[4]]
        
         }
    }
    
    for (or in 1:run1)
        loo1[,or]<-loo1[,or][order(loo1.k[,or])]
    
    loo1_mean<-apply(loo1,1,mean)
    loo1_sd<-apply(loo1,1,sd)
    names(loo1_mean)<-row.names(train_set)
    names(loo1_sd)<-row.names(train_set)
    row.names(loo1)<-row.names(train_set)
    colnames(loo1)<-paste("run",1:run1)
    dif<-loo1_mean-train_env                             
    max_er.c<-max(abs(dif))                                   
    mean_er.c<-mean(dif)                                 
    RMSEP<-sqrt((sum((dif)^2))/dim(train_set)[1])
    R2.c<-summary(lm(loo1_mean~train_env))$r.squared
    

error<-list(R2.c,mean_er.c,max_er.c,RMSEP,loo1_mean,loo1_sd,loo1)
error

}  # end of 10cross


############################     boot train_set  ###################################

boot1.wa<-function(train_set,train_env,boot=run,...)
    {    
    loo3<-matrix(NA,nrow=(dim(train_set)[1]*boot),ncol=3)
    loo4<-NA
    nam1<-row.names(train_set)
    dim_ts1<-dim(train_set)[1]
    dim_ts2<-dim(train_set)[1]-1    
    sample_n<-c(1:dim_ts1)  
    jb<-0
    result_inf<-rep(NA,dim_ts1)
    result_sd<-rep(NA,dim_ts1)
    set.seed(seed)
    for (n in seq_len(boot))
            {
            k<-sample((1:dim_ts1),replace=TRUE)
            sample_n1<-sample_n[sample_n%in%k=="FALSE"]
            l<-length(sample_n1)
            jb1<-c((jb+1):(jb+l))
            jb<-jb+l            
            train_set.b<-train_set[k[order(k)],]
            train_env.b<-train_env[k[order(k)]]            
            test_set.c<-train_set[sample_n1,]
            train_set.b<-train_set.b[,colSums(train_set.b)%in%0=="FALSE"]
            loo3[jb1,1]<-wa1(train_set.b,train_env.b,test_set.c)[[4]]
            loo3[jb1,2]<-train_env[sample_n1,]           
            loo4[jb1]<-sample_n1
            }
    loo3<-as.data.frame(loo3)
    loo3[1:length(loo4),3]<-loo4
    loo3<-loo3[1:length(loo4),]
    for (k in seq_len(dim_ts1))
            {
            result_inf[k]<-mean(loo3[1:jb,1][loo3[1:jb,3]==k],na.rm=TRUE)
            if (any(loo3[,3]==k))
                result_sd[k]<-sd(loo3[1:jb,1][loo3[1:jb,3]==k],na.rm=TRUE)
            }
    names(result_inf)<-row.names(train_set)
    names(result_sd)<-rownames(train_set)
    s1<-result_sd
    s2<-sqrt(sum((result_inf-train_env[,1])^2)/length(train_env[,1]))
    pred.error<-result_inf-train_env[,1]
    ms_s1 <- sqrt(mean(s1^2,na.rm=TRUE))
    ms_s2 <- sqrt(mean(pred.error^2,na.rm=TRUE))
    s_rmsep <- sqrt(s1^2 + ms_s2^2)
    ms_rmsep <- sqrt(ms_s1^2 + ms_s2^2)
    mean_error<-result_inf-train_env
    mean_error1<-mean(loo3[,1]-loo3[,2],na.rm=TRUE)
    max_error<-max(abs(result_inf-train_env),na.rm=TRUE)
    R2.c<-summary(lm(result_inf~train_env))$r.squared
    result<-list(R2.c,mean_error1,max_error,ms_rmsep,result_inf,mean_error,s1,s2,ms_s1,ms_s2,s_rmsep,result_sd)
}

boot2.wa<-function(train_set,train_env,test_set,boot=run,...)
    {    
    loo3<-matrix(NA,nrow=(dim(train_set)[1]*boot),ncol=3)
    loo4<-NA
    loo_test<-matrix(NA,nrow=dim(test_set)[1],ncol=boot)
    nam1<-row.names(train_set)
    dim_ts1<-dim(train_set)[1]
    dim_ts2<-dim(train_set)[1]-1    
    sample_n<-c(1:dim_ts1)  
    jb<-0
    result_inf<-rep(NA,dim_ts1)
    result_sd<-rep(NA,dim_ts1)
    set.seed(seed)
    for (n in seq_len(boot))
            {
            k<-sample((1:dim_ts1),replace=TRUE)
            sample_n1<-sample_n[sample_n%in%k=="FALSE"]
            l<-length(sample_n1)
            jb1<-c((jb+1):(jb+l))
            jb<-jb+l            
            train_set.b<-train_set[k[order(k)],]
            train_env.b<-train_env[k[order(k)]]            
            test_set.c<-train_set[sample_n1,]
            train_set.b<-train_set.b[,colSums(train_set.b)%in%0=="FALSE"]
            mod<-wa1(train_set.b,train_env.b,test_set.c)
            loo3[jb1,1]<-mod[[4]]
            opt<-mod[[8]]
            k<-colnames(test_set) %in% names(opt)
            k<-colnames(test_set)[k]
            test_set.opt<-test_set[,k]                                                       #Skalierung auf 100%???
            opt_d<-opt[k]
            if(dim(test_set)[1]>1)
                loo_test[,n]<-colSums(opt_d*t(test_set.opt),na.rm=TRUE)/rowSums(test_set.opt)      
            if(dim(test_set)[1]==1)
                loo_test[,n]<-sum(opt_d*t(test_set.opt),na.rm=TRUE)/sum(test_set.opt) 
            loo3[jb1,2]<-train_env[sample_n1,]           
            loo4[jb1]<-sample_n1
            }
    loo3<-as.data.frame(loo3)
    loo3[1:length(loo4),3]<-loo4
    loo3<-loo3[1:length(loo4),]
    for (k in seq_len(dim_ts1))
            {
            result_inf[k]<-mean(loo3[1:jb,1][loo3[1:jb,3]==k],,na.rm=TRUE)
            if (any(loo3[,3]==k))
                result_sd[k]<-sd(loo3[1:jb,1][loo3[1:jb,3]==k],na.rm=TRUE)
            }
s1<-result_sd
s2<-sqrt(sum((result_inf-train_env[,1])^2)/length(train_env[,1]))
pred.error<-result_inf-train_env[,1]
ms_s1 <- sqrt(mean(s1^2,na.rm=TRUE))
ms_s2 <- sqrt(mean(pred.error^2,na.rm=TRUE))
s_rmsep <- sqrt(s1^2 + ms_s2^2)
ms_rmsep <- sqrt(ms_s1^2 + ms_s2^2)
mean_error<-result_inf-train_env
mean_error1<-mean(loo3[,1]-loo3[,2],na.rm=TRUE)
max_error<-max(abs(result_inf-train_env),na.rm=TRUE)

R2.c<-summary(lm(result_inf~train_env))$r.squared
loo_test.mean<-apply(loo_test,1,function(x) mean(x,na.rm=TRUE))
loo_test.sd<-apply(loo_test,1,function(x) sd(x,na.rm=TRUE))
row.names(loo_test)<-row.names(test_set)
names(loo_test.mean)<-row.names(test_set)
names(loo_test.sd)<-row.names(test_set)
names(result_sd)<-row.names(train_set)
names(result_inf)<-row.names(train_set)
result<-list(R2.c,mean_error1,max_error,ms_rmsep,result_inf,mean_error,s1,s2,ms_s1,ms_s2,s_rmsep,result_sd,loo_test.mean,loo_test.sd,loo_test)
result
}
#####################################   end of functions       #############################################


############################################################################################################


#####################################    MAIN CODE             #############################################
if (data_l==3)
        run_wa<-wa1(train_set,train_env,test_set)
if (data_l==2)
        run_wa<-wa1(train_set,train_env)
run_wa1<-wa1(train_set,train_env)

if (val=="loo")
    val1<-"Leave-one-out"
if (val=="n.cross") 
    val1<-"n-fold-cross validation"
if (val=="boot")
    val1="bootstrap"

if (out=="TRUE")
 {
    cat("",fill=TRUE)
    cat("",fill=TRUE)
    cat("                   transfer function",fill=TRUE)
    cat("",fill=TRUE)
    cat("",fill=TRUE)
    cat("type        = weighted averaging",fill=TRUE)
    cat("method      =",desh.meth,fill=TRUE)
    cat("n samples   =",dim(train_set)[1],fill=TRUE)
    cat("n species   =",dim(train_set)[2],fill=TRUE)
    cat("val.-method =",val1,fill=TRUE) 
    if (val=="n.cross")
        {             
            cat("       seed =",seed,fill=TRUE)
            cat("        run =",run,fill=TRUE)
            cat("     n.fold =",nfold,fill=TRUE)
        }
    if (val=="boot")
         {             
            cat("       seed =",seed,fill=TRUE)
            cat("        run =",run,fill=TRUE)
        }    
    cat("",fill=TRUE)
   cat("R2    = ",run_wa1[[5]],"          Mean-error = ",round(run_wa1[[2]],4),fill=TRUE)
   cat("RMSE  = ",run_wa1[[6]],"          Max-error  = ",round(run_wa1[[3]],4),fill=TRUE)
   if (val!="none")
            {
            cat("",fill=TRUE)
            cat("                     * please wait *",fill=TRUE) 
             }
 flush.console()
 }

error<-list(NA,NA,NA,NA,NA,NA,NA,NA)

if (val=="loo")
    {   if (data_l==2)
            error<-loo.wa(train_set,train_env)
        if (data_l==3)
            error<-loo_pred.wa(train_set,train_env,test_set)
    }    
if (val=="n.cross")  
            error<-tencross.wa(train_set,train_env)
if (val=="boot")
     {   if (data_l==2)
            error<-boot1.wa(train_set,train_env)
         if (data_l==3)
            error<-boot2.wa(train_set,train_env,test_set)
     }

if (desh.meth=="inverse")
error_m<-matrix(nrow=1,ncol=8,dimnames=list("wa-inv.desh.",c("RMSE","R2","Ave_Bias","Max_Bias","X_R2","X_Ave_Bias","X_Max_Bias","RMSEP")))
if (desh.meth=="class")
error_m<-matrix(nrow=1,ncol=8,dimnames=list("wa-class.desh.",c("RMSE","R2","Ave_Bias","Max_Bias","X_R2","X_Ave_Bias","X_Max_Bias","RMSEP")))

error_m[1,1]<-run_wa1[[6]]
error_m[1,2]<-run_wa1[[5]]
error_m[1,3]<-run_wa1[[2]]
error_m[1,4]<-run_wa1[[3]]
error_m[1,5]<-error[[1]]
error_m[1,6]<-error[[2]]
error_m[1,7]<-error[[3]]
error_m[1,8]<-error[[4]]

if (data_l==2)
{
    if (val=="none")
        results<-list(spec_n_train,n2_train,run_wa1[[1]],run_wa1[[7]],error_m)
    
    
    
    if (val=="loo")
        {
        results<-list(spec_n_train,n2_train,run_wa1[[1]],run_wa1[[7]],error_m,error[[5]])
        names(results)[[6]]<-"inferred train.set (loo)"
         }
    if (val=="n.cross")
        {
        results<-list(spec_n_train,n2_train,run_wa1[[1]],run_wa1[[7]],error_m,error[[5]],error[[6]],error[[7]])
        names(results)[[6]]<-"mean(inferred train.set) (n.cross)"
        names(results)[[7]]<-"sd(inferred train.set) (n.cross)"
        names(results)[[8]]<-"inferred train.set (n.cross)"
        }
    
    if (val=="boot")
        {
        results<-list(spec_n_train,n2_train,run_wa1[[1]],run_wa1[[7]],error_m,error[[5]],error[[11]],error[[9]],error[[10]])
        names(results)[[6]]<-"mean(inferred train.set) (boot)"
        names(results)[[7]]<-"sd(inferred train.set) (boot)"
        names(results)[[8]]<-"s1"
        names(results)[[9]]<-"s2"
        }
    names(results)[[1]]<-"species in train.set"
    names(results)[[2]]<-"N2 train.set"
    names(results)[[3]]<-"species.optima"
    names(results)[[4]]<-"inferred train.set"
    names(results)[[5]]<-"performance"
}

if (data_l==3)
{
    if (val=="none")
        results<-list(spec_n_train,n2_train,run_wa1[[1]],run_wa1[[7]],error_m,spec_n_test,n_test.train,n2_test,run_wa[[4]])
    if (val=="loo")
        results<-list(spec_n_train,n2_train,run_wa1[[1]],run_wa1[[7]],error_m,error[[5]],spec_n_test,n_test.train,n2_test,run_wa[[4]],error[[6]],error[[7]],error[[8]])
    if (val=="n.cross")
        results<-list(spec_n_train,n2_train,run_wa1[[1]],run_wa1[[7]],error_m,error[[5]],error[[6]],spec_n_test,n_test.train,n2_test,run_wa[[4]])
    if (val=="boot")
        results<-list(spec_n_train,n2_train,run_wa1[[1]],run_wa1[[7]],error_m,error[[5]],error[[12]],error[[9]],error[[10]],spec_n_test,n_test.train,n2_test,run_wa[[4]],error[[13]],error[[14]],error[[15]])
    
    names(results)[[1]]<-"species in train.set"
    names(results)[[2]]<-"N2 train.set"
    names(results)[[3]]<-"species.optima"
    names(results)[[4]]<-"inferred train.set"
    names(results)[[5]]<-"performance"
    if (val=="none")    
        {
        names(results)[[6]]<-"species in core.samples"
        names(results)[[7]]<-"n species core.samples in train.set"
        names(results)[[8]]<-"N2 in core.samples"
        names(results)[[9]]<-"reconstruction_core.samples"
        }
    if (val=="loo")
        {
        names(results)[[6]]<-"inferred train.set.val"
        names(results)[[7]]<-"species in core.samples"
        names(results)[[8]]<-"n species core.samples in train.set"
        names(results)[[9]]<-"N2 in core.samples"
        names(results)[[10]]<-"reconstruction_core.samples"
        names(results)[[11]]<-"mean(reconstruction_core.samples).val"
        names(results)[[12]]<-"sd(reconstruction_core.samples).val"
        names(results)[[13]]<-"reconstruction_core.samples.val"
        }
    if (val=="n.cross")
        {
        names(results)[[6]]<-"inferred train.set.val"
        names(results)[[7]]<-"sd(reconstruction_core.samples).val"
        names(results)[[8]]<-"species in core.samples"
        names(results)[[9]]<-"n species core.samples in train.set"
        names(results)[[10]]<-"N2 in core.samples"
        names(results)[[11]]<-"reconstruction_core.samples"        
        }

    if (val=="boot")
        {names(results)[[6]]<-"mean(inferred train.set).val"
        names(results)[[7]]<-"sd(inferred train.set).val"
        names(results)[[8]]<-"s1"
        names(results)[[9]]<-"s2"
        names(results)[[10]]<-"species in core.samples"
        names(results)[[11]]<-"n species core.samples in train.set"
        names(results)[[12]]<-"N2 in core.samples"
        names(results)[[13]]<-"reconstruction_core.samples"
        names(results)[[14]]<-"mean(reconstruction_core.samples).val"
        names(results)[[15]]<-"sd(reconstruction_core.samples).val"
        names(results)[[16]]<-"reconstruction_core.samples.val"
        }
}

if (out=="TRUE" && val!="none")
    {
        
      cat("",fill=TRUE)
   
         if (val=="boot")
            {cat("s1       = ",error[[9]],  "          s2               = ",round(error[[10]],4),fill=TRUE)}
            cat("R2.c     = ",error_m[1,5],"          Mean-error.c     = ",round(error_m[1,6],4),fill=TRUE)
            cat("RMSEP    = ",error_m[1,8],"          Max-error.c      = ",round(error_m[1,7],4),fill=TRUE)
        cat("",fill=TRUE)
        cat("",fill=TRUE)
    }    



##################################################################                PLOTS                   ###################################

if(d.plot==TRUE)
    {
    if (val=="none")
        {   
        par(mfrow=c(1,2))
        x1<-range(train_env)[1]-range(train_env)[1]*0.5
        x2<-range(train_env)[2]+range(train_env)[2]*0.2
        plot(run_wa1[[7]]~train_env,xlab="observed",ylab="inferred (train_set)",col="gray",pch=19,xlim=c(x1,x2),ylim=c(x1,x2))
        abline(0,1,col="red",lty=2)
        panel.smooth(train_env,run_wa1[[7]],col.smooth="blue",lwd=2)
        dif.plot<-(run_wa1[[7]]-train_env)
        plot(dif.plot~train_env,col="gray",pch=19,xlim=c(x1,x2),xlab="observed",ylab="residuals (inferred-measured)")
        abline(0,0,col="red",lty=2)
        panel.smooth(train_env,dif.plot,col.smooth="blue",lwd=2)
        }
    if (val!="none")
        {   
        par(mfrow=c(2,2))
        x1<-range(train_env)[1]-range(train_env)[1]*0.5
        x2<-range(train_env)[2]+range(train_env)[2]*0.2
        plot(run_wa1[[7]]~train_env,xlab="observed",ylab="inferred (train_set)",col="gray",pch=19,xlim=c(x1,x2),ylim=c(x1,x2))
        abline(0,1,col="red",lty=2)
        panel.smooth(train_env,run_wa1[[7]],col.smooth="blue",lwd=2)
        dif.plot<-(run_wa1[[7]]-train_env)
        plot(dif.plot~train_env,col="gray",pch=19,xlim=c(x1,x2),xlab="observed",ylab="residuals (inferred-measured)")
        abline(0,0,col="red",lty=2)
        panel.smooth(train_env,dif.plot,col.smooth="blue",lwd=2)
        plot(error[[5]]~train_env,xlab="observed",ylab=paste("inferred  (",val1,")"),col="gray",pch=19,xlim=c(x1,x2),ylim=c(x1,x2))
        abline(0,1,col="red",lty=2)
        panel.smooth(train_env,error[[5]],col.smooth="blue",lwd=2)
        dif.plot<-(error[[5]]-train_env)
        plot(dif.plot~train_env,col="gray",pch=19,xlim=c(x1,x2),xlab="observed",ylab="residuals (cross.validation-measured)")
        abline(0,0,col="red",lty=2)
        panel.smooth(train_env,dif.plot,col.smooth="blue",lwd=2)
        }
    }
#################################################################################################################################################

invisible(results)

}

