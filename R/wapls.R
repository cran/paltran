

wapls <-function(...,comp=4,d.plot=TRUE,plot.comp="RMSEP",env.trans=FALSE,spec.trans=FALSE,diagno=TRUE,seed=1,run=10,val=c("none","10-cross","loo","boot"),scale=FALSE,out=TRUE,drop.non.sig=FALSE,min.occ=1)
{

    n_comp<-comp
    data<-list(...)
    train_set<-as.matrix(data[[1]])
    train_env<-as.matrix(data[[2]])
    train_set<-train_set[,colSums(train_set)!=0]
    if (missing(val)) 
        {
        val <- "none"
        val1<- "none"
        }
    val<- match.arg(val)

    scores_wapls<-NA
    inf_train_wapls<-NA 
    inf_test_wapls<-NA
    error<-NA

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
            for (tax in 1:col)
                {
                    env<-train_env[,1]
                    if (sum(ifelse(train_set[,tax]>0,1,0))>5)
                        {
                            fit.gam<-gam(train_set[,tax]~s(env))    
                            k<-summary(fit.gam)[[8]]
                            if (k<0.05)
                                sig.ts[tax]<-1 else sig.ts[tax]<-0
                        } else    sig.ts[tax]<-0 
                }
             train_set<-train_set[,sig.ts==1]
            train_set1<-train_set[rowSums(train_set)!=0,]
            train_env<-as.matrix(train_env[rowSums(train_set)!=0,])
            train_set<-train_set1
        }
  
   
 
 
    dat<-length(data)
    if (dat==3)
        test_set<-data[[3]]
    if (dat<=2)
        test_set<-NA
  
    col1<-dim(train_set)[2]
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
    n_test.train<-NA
if (diagno==TRUE)
    {
        
        library(vegan)  
        spec_n_train<-specnumber(train_set)
        n2_train<-1/rowSums((train_set/100)^2)
        if (length(data)==3)
            if(dim(test_set)[1]>=2)
            {
                spec_n_test<-specnumber(test_set)
                n2_test<-1/rowSums((test_set/100)^2)
                k<-colnames(test_set) %in% colnames(train_set)
                k<-colnames(test_set)[k]
                test_set<-test_set[,k]                 
                n_test.train<-specnumber(test_set)
            }
    }

sort0<-function(x)
    {
        x<-as.data.frame(x)
        x<-x[,order(-colSums(x))]
        repeat
        { 
            if (sum(x[,dim(x)[2]])==0)
            x<-x[,-dim(x)[2]] 
            if (sum(x[,dim(x)[2]])!=0) break
        } 
    x
    }
train_set<-sort0(train_set)
if (scale==TRUE)
    {
    if (length(data)==3)
        {test_set<-data[[3]]
        if (dim(test_set)[1]>1)
            test_set<-as.data.frame(t(apply(test_set,1,function(x)x*100/sum(x,na.rm=TRUE)))) 
        }

    train_set<-as.data.frame(t(apply(train_set,1,function(x)x*100/sum(x,na.rm=TRUE))))
    }            

###########################  declaration of variables #########################
RMSE<-NA
mean_error<-NA
max_error<-NA
R2<-NA
u<-NA
col<-dim(train_set)[2]
row<-dim(train_set)[1]
numb.c<-c(1:n_comp)
env_inf_train<-matrix(ncol=n_comp,nrow=row)
env_infp<-matrix(0,ncol=n_comp,nrow=row1)
nam<-c("comp1","comp2","comp3","comp4","comp5","comp6","comp7","comp8","comp9","comp10","comp11","comp12","comp13","comp14")
nam1<-nam[1:n_comp]
error<-matrix(0,nrow=n_comp,ncol=10,dimnames=list(nam1,c("RMSE","R2","Ave_Bias","Max_Bias","X_R2","X_Ave_Bias","X_Max_Bias","RMSEP","n_comp","MW")))
fac<-c(1:n_comp)
env_infd.cross<-matrix(nrow=row,ncol=n_comp)
#############

wapls1<-function(...,n_comp=comp)

{    
    data<-list(...)
    train_set<-data[[1]]
    train_env<-as.matrix(data[[2]])
    dat<-length(data)
    row1<-1
    u<-NA
    z<-NA
    env_inf<-NA
    r.n1.pre<-NA
    r.n1.pre<-NA
    score<-NA
    RMSE<-NA
    mean_error<-NA
    max_error<-NA
    R2<-NA
    if (dat==3)
    {        
        test_set<-data[[3]] 
        if (dim(test_set)[1]>1)
            {sum2<-rowSums(test_set)  
             row1<-dim(test_set)[1]    
            }
        if (dim(test_set)[1]==1)
            {sum2<-sum(test_set)  
             row1<-1    
            }   
        r.n1.pre<-matrix(ncol=n_comp,nrow=row1)
        env_inf<-matrix(0,ncol=n_comp,nrow=row1)
        r.n1.pre<-matrix(ncol=n_comp,nrow=row1)   
    }
    sum_sd<-sum(train_env*rowSums(train_set),na.rm=TRUE)/sum(rowSums(train_set,na.rm=TRUE))
    train_env.c<-train_env-sum_sd
    r<-train_env.c 
    train_env1.c<-as.matrix(train_env.c)
    col<-ncol(train_set)
    row<-nrow(train_set)
    r.n1<-matrix(ncol=n_comp,nrow=row)
    env_inf_train<-matrix(0,ncol=n_comp,nrow=row)
    train_set.rs<-rowSums(train_set,na.rm=TRUE)
    updated.opt<-matrix(0,ncol=n_comp,nrow=col)
    row.names(updated.opt)<-names(train_set)
    divid<-sum(train_set.rs)
    r.n_all<-matrix(NA,ncol=n_comp,nrow=row)
    for (run in 1:n_comp)
    {
        for(i in seq_len(col))
            u[i]<-sum(r*train_set[,i],na.rm=TRUE)/sum(train_set[,i],na.rm=TRUE)       
        r.n<-as.vector(colSums(u*t(train_set),na.rm=TRUE)/train_set.rs)
        updated.opt[,run]<-u+sum_sd
        if (run>1)
            {
                v<-sum(((train_set.rs*r.n)*r.n1[,(run-1)]))/sum(train_set.rs)
                r.n<-r.n-(v*r.n1[,run-1])
            }
    
        r.n<-as.vector(r.n)
        r.n_all[,run]<-r.n
        z<-sum(r.n*train_set.rs)/divid
            s2<-sum(((r.n-z)^2)*train_set.rs)/divid
                r.n1[,run]<-(r.n-z)/sqrt(s2)

        score<-r.n1[,1:run]
        weights1<-train_set.rs/divid
        fit.lm<-lm(train_env1.c~score,weights=(weights1))

        env_inf_train.c<-fitted(fit.lm)
        env_inf_train[,run]<-env_inf_train.c+sum_sd
        r<-as.vector((-1)*(train_env1.c-env_inf_train.c))

        dif<-(env_inf_train[,run]-train_env)
        y<-as.matrix(train_env)
        RMSE[run]<-sqrt(sum((dif)^2)/nrow(train_set))
        mean_error[run]<-mean((dif))
        max_error[run]<-max(abs(dif))
        R2[run]<-summary(lm(y~env_inf_train[,run]))$r.squared
# prediction for test_set
        if (dat==3)
            {
            k<-colnames(test_set) %in% colnames(train_set)       
            k<-colnames(test_set)[k]
            test_set.n<-test_set[,k]
            train_set.n<-train_set[k]
            u.v<-as.vector(u)
            names(u.v)<-names(train_set)
            u.v<-u.v[k]
            r.n.pre<-NA
                if (dim(test_set)[1]>1)
                    r.n.pre<-colSums(u.v*t(test_set.n),na.rm=TRUE)/rowSums(test_set.n,na.rm=TRUE)
                if (dim(test_set)[1]==1)
                    r.n.pre<-sum(u.v*t(test_set.n),na.rm=TRUE)/sum(test_set.n,na.rm=TRUE)
                
                if (run>1)
                    r.n.pre<-r.n.pre-(v*r.n1.pre[,run-1])
                r.n1.pre[,run]<-(r.n.pre-z)/sqrt(s2)
                for (k in 1:run)
                    env_inf[,run]<-env_inf[,run]+fit.lm$coefficients[k+1]*r.n1.pre[,k]
                env_inf[,run]<-env_inf[,run]+fit.lm$coefficients[1]
                env_inf[,run]<-env_inf[,run]+sum_sd
            }      
    }
    result<-list(r,R2,RMSE,mean_error,max_error,env_inf,r.n_all,env_inf_train,r.n1,u,r.n1,updated.opt)
  result  
}   #end of wapls1

##############################         loo                                           ###########################


loo.wapls<-function(train_set,train_env,n_comp=comp)
{
    dim_row<-nrow(train_set)
    loo1<-matrix(nrow=dim_row,ncol=n_comp)
    error1<-matrix(nrow=n_comp,ncol=4)
    train_env<-as.matrix(train_env)
    for (kn in seq_len(dim_row))
        {
        train_set.c<-train_set[-kn,]
        train_set.c<-train_set.c[,order(-colSums(train_set.c))]
        dim_ts<-ncol(train_set.c)
        train_set.c<-train_set.c[,colSums(train_set.c)!=0]                                 #######    !!!!! in C2 nicht eingeschaltet!!!!!!
        loo1[kn,]<-wapls1(train_set.c,train_env[-kn,],train_set[kn,])[[6]]
        }
    dif<-(loo1-train_env[,1])
    for (k in seq_len(n_comp))
    {
    error1[k,3]<-max(abs(dif[,k]))
    error1[k,2]<-mean(dif[,k])
    error1[k,4]<-sqrt(sum((dif[,k])^2)/dim_row)
    error1[k,1]<-summary(lm(loo1[,k]~train_env[,1]))$r.squared
    }
    result<-list(error1[,1],error1[,2],error1[,3],error1[,4],loo1)
    
result
}

loo2.wapls<-function(train_set,train_env,test_set,n_comp=comp)
{
    dim_row<-nrow(train_set)
    loo1<-matrix(nrow=dim_row,ncol=n_comp)
    dim_test<-dim(test_set)[1]
    loo_pred<-matrix(nrow=dim_row*dim_test,ncol=n_comp+1)
    error1<-matrix(nrow=n_comp,ncol=4)
    train_env<-as.matrix(train_env)
    result_mean.test<-matrix(nrow=dim_test,ncol=n_comp)
    result_sd.test<-matrix(nrow=dim_test,ncol=n_comp)
    for (kn in seq_len(dim_row))
        {
        train_set.c<-train_set[-kn,]
        train_set.c<-train_set.c[,order(-colSums(train_set.c))]
        dim_ts<-ncol(train_set.c)
        train_set.c<-train_set.c[,colSums(train_set.c)!=0]                                 #######    !!!!! in C2 nicht eingeschaltet!!!!!!
        loo1[kn,]<-wapls1(train_set.c,train_env[-kn,],train_set[kn,])[[6]]
        s1<-(kn-1)*(dim_test)+1
        s2<-kn*dim_test
        loo_pred[s1:s2,1:n_comp]<-wapls1(train_set.c,train_env[-kn,],test_set)[[6]]
        loo_pred[s1:s2,n_comp+1]<-1:dim_test
        }
    dif<-(loo1-train_env[,1])
    for (k in seq_len(n_comp))
    {
    error1[k,3]<-max(abs(dif[,k]))
    error1[k,2]<-mean(dif[,k])
    error1[k,4]<-sqrt(sum((dif[,k])^2)/dim_row)
    error1[k,1]<-summary(lm(loo1[,k]~train_env[,1]))$r.squared
    }
    
    
        for (k in seq_len(dim_test))
            {
            result_mean.test[k,]<-apply(loo_pred[loo_pred[,n_comp+1]==k,1:n_comp],2,mean)
            result_sd.test[k,]<-apply(loo_pred[loo_pred[,n_comp+1]==k,1:n_comp],2,sd)
            }
    
    
    result<-list(error1[,1],error1[,2],error1[,3],error1[,4],loo1,result_mean.test,result_sd.test)
    
result
}

##############################     cross_vall                                          ###########################

tencross.wapls<-function(train_set,train_env,run1=run,n_comp=comp)
{  
   max_er.c<-NA
   mean_er.c<-NA
   RMSEP<-NA
   c.cross<-NA
   R2.c<-NA
   dim_row<-nrow(train_set)
   loo1<-matrix(nrow=dim_row*run1,ncol=n_comp+1)
   loo2<-vector("list",length=n_comp)
   result_mean.train<-matrix(nrow=dim_row,ncol=n_comp)
   result_sd.train<-matrix(nrow=dim_row,ncol=n_comp)
   error1<-matrix(NA,nrow=n_comp,ncol=4)
   inf_train.cross<-matrix(nrow=dim_row,ncol=n_comp)
   cross_number<-round(dim_row/10,0)
   error_max<-matrix(ncol=n_comp,nrow=run1)
   error_mean<-matrix(ncol=n_comp,nrow=run1)
   error_RMSEP<-matrix(ncol=n_comp,nrow=run1)
   error_R2<-matrix(ncol=n_comp,nrow=run1)
   set.seed(seed)
    for (i in seq_len(10))
        c.cross[i]<-(i-1)*cross_number
    c.cross[11]<-max(length(train_env))
   for (r in seq_len(run1))
    {      
        k<-sample((1:dim(train_set)[1]))
        train_set.b<-train_set[k,]
        train_env.b<-train_env[k,]
        for (ic in seq_len(10))
        {
            c1<-c.cross[ic]+1
            c2<-c.cross[ic+1]
            c1.l<-c1+((r-1)*dim_row)
            c2.l<-c2+((r-1)*dim_row)
            l_train<-nrow(train_env)
            train_set.c<-train_set.b[-(c1:c2),]
            train_env.c<-train_env.b[-(c1:c2)]
            test_set.c<-train_set.b[(c1:c2),]       
            train_set.c<-train_set.c[,order(-colSums(train_set.c))]
            train_set.c<-train_set.c[,colSums(train_set.c)!=0]
            loo1[c1.l:c2.l,1:n_comp]<-wapls1(train_set.c,train_env.c,test_set.c)[[6]]        
            loo1[c1.l:c2.l,n_comp+1]<-k[c1:c2]
        }
    } 
   for (k in seq_len(dim_row))
        {
            result_mean.train[k,]<-apply(loo1[loo1[,n_comp+1]==k,1:n_comp],2,mean)
            result_sd.train[k,]<-apply(loo1[loo1[,n_comp+1]==k,1:n_comp],2,sd)
        }
   for (i in seq_len(n_comp))
        {
            dif<-result_mean.train[,i]-train_env
            error1[i,1]<-summary(lm(result_mean.train[,i]~train_env))$r.squared
            error1[i,2]<-mean(dif)
            error1[i,3]<-max(abs(dif))
            error1[i,4]<-sqrt(mean(dif^2))
         }     
    error<-list(error1[,1],error1[,2],error1[,3],error1[,4],result_mean.train,result_sd.train)
   
} 
#try<-tencross.wapls(train_set,train_env,run1=10)


###############################################    bootstrap  ##########################

boot1.wapls<-function(train_set,train_env,boot=run,...)
{    
 #boot=100
    loo3<-matrix(NA,nrow=(dim(train_set)[1]*boot),ncol=n_comp+2)
    #loo5<-matrix(NA,nrow=(dim(test_set)[1]*boot),ncol=n_comp+1)
    loo4<-NA
    loo7<-NA
    nam1<-row.names(train_set)
    loo3.res<-matrix(ncol=3,nrow=dim(train_set)[1],dimnames=list(nam1,c("under b-error","mean","upper b-error")))   
    dim_ts1<-dim(train_set)[1]
    dim_ts2<-dim(train_set)[1]-1    
    result_inf<-matrix(nrow=dim_ts1,ncol=n_comp)    
    result_sd<-matrix(nrow=dim_ts1,ncol=n_comp) 
    train_env<-as.matrix(train_env)
    spec_n<-c(1:dim_ts1)  
    jb<-0
    #ntest<-dim(test_set)[1]
    
    set.seed(seed)
    for (n in seq_len(boot))
            {
            k<-sample((1:dim_ts1),replace=TRUE)
            spec_n1<-spec_n[spec_n%in%k=="FALSE"]
            l<-length(spec_n1)
            jb1<-c((jb+1):(jb+l))
            jb<-jb+l            
            train_set.b<-train_set[k[order(k)],]
            train_env.b<-train_env[k[order(k)]]            
            test_set1<-train_set[spec_n1,]
            train_set.b<-train_set.b[,order(-colSums(train_set.b))]      
            train_set.b<-train_set.b[,colSums(train_set.b)%in%0=="FALSE"]
            loo3[jb1,1:n_comp]<-wapls1(train_set.b,train_env.b,test_set1)[[6]]
            loo3[jb1,n_comp+1]<-train_env[spec_n1,]           
            loo4[jb1]<-spec_n1
            
            }
    loo3<-as.data.frame(loo3)
    loo3[1:length(loo4),n_comp+2]<-loo4
    nam_ds<-levels(as.factor(loo3[,n_comp+2]))
    for (n in 1:comp)
        for (k in seq_len(dim_ts1))
            {
            result_inf[k,n]<-mean(loo3[1:jb,n][loo3[1:jb,n_comp+2]==k],,na.rm=TRUE) 
            if (any(loo3[,n_comp+2]==k,na.rm=TRUE))
                result_sd[k,n]<-sd(loo3[1:jb,n][loo3[1:jb,n_comp+2]==k],na.rm=TRUE)     
            }
   

s1<-result_sd
s2<-sqrt(colSums((result_inf-train_env[,1])^2,na.rm=TRUE)/length(train_env[,1]))
pred.error<-result_inf-train_env[,1]
ms_s1 <- sqrt(colMeans(s1^2,na.rm=TRUE))
ms_s2 <- sqrt(colMeans(pred.error^2,na.rm=TRUE))
s_rmsep <- sqrt(s1^2 + ms_s2^2)
ms_rmsep <- sqrt(ms_s1^2 + ms_s2^2)
R2.c<-rep(NA,n_comp)
for (r in 1:n_comp)
    R2.c[r]<-summary(lm(result_inf[,r]~train_env[,1]))$r.square

mean_error<-apply((result_inf-train_env[,1]),2,function(x) mean(x,na.rm=TRUE))
max_error<-apply((result_inf-train_env[,1]),2,function(x) max(abs(x),na.rm=TRUE))
result<-list(R2.c,mean_error,max_error,ms_rmsep,result_inf,ms_s1,ms_s2,result_sd)
result
}

boot2.wapls<-function(train_set,train_env,test_set,boot=run,...)
{    
    #boot=100
    loo3<-matrix(NA,nrow=(dim(train_set)[1]*boot),ncol=n_comp+2)
    loo5<-matrix(NA,nrow=(dim(test_set)[1]*boot),ncol=n_comp+1)
    loo4<-NA
    loo7<-NA
    nam1<-row.names(train_set)
    loo3.res<-matrix(ncol=3,nrow=dim(train_set)[1],dimnames=list(nam1,c("under b-error","mean","upper b-error")))   
    dim_ts1<-dim(train_set)[1]
    dim_ts2<-dim(train_set)[1]-1    
    result_inf<-matrix(nrow=dim_ts1,ncol=n_comp)    
    result_sd<-matrix(nrow=dim_ts1,ncol=n_comp) 
    train_env<-as.matrix(train_env)
    spec_n<-c(1:dim_ts1)  
    jb<-0
    ntest<-dim(test_set)[1]
    result_inf_ts<-matrix(nrow=ntest,ncol=n_comp) 
    result_sd_ts<-matrix(nrow=ntest,ncol=n_comp) 
    set.seed(seed)
    for (n in seq_len(boot))
            {
            k<-sample((1:dim_ts1),replace=TRUE)
            spec_n1<-spec_n[spec_n%in%k=="FALSE"]
            l<-length(spec_n1)
            jb1<-c((jb+1):(jb+l))
            jb<-jb+l            
            jb2<-c((((n-1)*ntest)+1):(n*ntest))
            train_set.b<-train_set[k[order(k)],]
            train_env.b<-train_env[k[order(k)]]            
            test_set1<-train_set[spec_n1,]
            train_set.b<-train_set.b[,order(-colSums(train_set.b))]      
            train_set.b<-train_set.b[,colSums(train_set.b)%in%0=="FALSE"]
            loo3[jb1,1:n_comp]<-wapls1(train_set.b,train_env.b,test_set1)[[6]]
            loo3[jb1,n_comp+1]<-train_env[spec_n1,]           
            loo4[jb1]<-spec_n1
            loo5[jb2,1:n_comp]<-wapls1(train_set.b,train_env.b,test_set)[[6]]
            loo7[jb2]<-c(1:length(row.names(test_set)))
            }
    loo3<-as.data.frame(loo3)
    loo3[1:length(loo4),n_comp+2]<-loo4
    nam_ds<-levels(as.factor(loo3[,n_comp+2]))
    loo5<-as.data.frame(loo5)
    loo5[1:length(loo7),n_comp+1]<-loo7
    nam_ts<-levels(as.factor(loo5[,n_comp+1]))
    for (n in 1:comp)
        for (k in seq_len(dim_ts1))
            {
            result_inf[k,n]<-mean(loo3[1:jb,n][loo3[1:jb,n_comp+2]==k],na.rm=TRUE)
            if (any(loo3[,n_comp+2]==k,na.rm=TRUE))
                result_sd[k,n]<-sd(loo3[1:jb,n][loo3[1:jb,n_comp+2]==k],na.rm=TRUE)
            }
    for (n in 1:comp)
        for (k in seq_len(ntest))
            {
            result_inf_ts[k,n]<-mean(loo5[1:jb,n][loo5[1:jb,n_comp+1]==k],,na.rm=TRUE)
            if (any(loo5[,n_comp+1]==k,na.rm=TRUE))
                result_sd_ts[k,n]<-sd(loo5[1:jb,n][loo5[1:jb,n_comp+1]==k],na.rm=TRUE)
            }

s1<-result_sd
s2<-sqrt(colSums((result_inf-train_env[,1])^2,na.rm=TRUE)/length(train_env[,1]))
pred.error<-result_inf-train_env[,1]
ms_s1 <- sqrt(colMeans(s1^2,na.rm=TRUE))
ms_s2 <- sqrt(colMeans(pred.error^2,na.rm=TRUE))
s_rmsep <- sqrt(s1^2 + ms_s2^2)
ms_rmsep <- sqrt(ms_s1^2 + ms_s2^2)
R2.c<-rep(NA,n_comp)
for (r in 1:n_comp)
    R2.c[r]<-summary(lm(result_inf[,r]~train_env[,1]))$r.square
mean_error<-apply((result_inf-train_env[,1]),2,function(x) mean(x,na.rm=TRUE))
max_error<-apply((result_inf-train_env[,1]),2,function(x) max(abs(x),na.rm=TRUE))
result<-list(result_sd_ts,s1,s2,result_inf_ts,result_inf,ms_s1,ms_s2,mean_error,max_error,s_rmsep,ms_rmsep,R2.c,result_sd)

result
}





###########################################################      main part    ##################################################


cval<-matrix(NA,ncol=4,nrow=comp)
cval_l<-NA
boot_inf_train<-NA
boot_recon<-NA
wapls_cross<-NA
boot_train<-NA



wapls_run<-wapls1(train_set,train_env)

if (dat==3)
    wapls_run2<-wapls1(train_set,train_env,test_set)
 
result1<-matrix(0,nrow=n_comp,ncol=10,dimnames=list(nam1,c("R2","Mean-error","Max_error","RMSE","R2.c","Mean-error.c","Max-error.c","RMSEP","n_comp","MW")))
        result1[,9]<-seq(1:n_comp)
        result1[,10]<-seq(1:n_comp)

result1[,1]<-wapls_run[[2]]
result1[,2]<-wapls_run[[4]]
result1[,3]<-wapls_run[[5]]
result1[,4]<-wapls_run[[3]]
scores_wapls<-round(wapls_run[[7]],6)
upd.opt<-round(wapls_run[[12]],3)
inf_train_wapls<-wapls_run[[8]]
inf_test_wapls<-wapls_run$reconstruction

if (val=="loo")
    val1<-"Leave-one-out"
if (val=="10-cross") 
    val1<-"10-fold-cross validation"
if (val=="boot")
    val1="bootstrap"

if (out=="TRUE")
{
 cat("",fill=TRUE)
 cat("",fill=TRUE)
 cat("                   transfer function",fill=TRUE)
 cat("",fill=TRUE)
 cat("",fill=TRUE)
 cat("type        = weighted averaging - partial least squares",fill=TRUE)
 cat("components  =",comp,fill=TRUE)
 cat("n samples   =",dim(train_set)[1],fill=TRUE)
 cat("n species   =",dim(train_set)[2],fill=TRUE)
 cat("val.-method =",val1,fill=TRUE) 
 if (val=="10-cross")
  {             
  cat("       seed =",seed,fill=TRUE)
  cat("        run =",run,fill=TRUE)
  }
  if (val=="boot")
   {             
   cat("       seed =",seed,fill=TRUE)
   cat("        run =",run,fill=TRUE)
   }    
   cat("",fill=TRUE)
   print(round(result1[,1:4],7))
   
   
   if (val!="none")
            {
            cat("",fill=TRUE)
            cat("                     * please wait *",fill=TRUE) 
             }
   flush.console()
}

if (dat==2)
{
    if (val=="none")
        results<-list(spec_n_train,n2_train,upd.opt,scores_wapls,inf_train_wapls,result1)

    
    
    if (val=="loo")
        {
        wapls_val<-loo.wapls(train_set,train_env)
        result1[,5]<-wapls_val[[1]]
        result1[,6]<-wapls_val[[2]]
        result1[,7]<-wapls_val[[3]]
        result1[,8]<-wapls_val[[4]]
        results<-list(spec_n_train,n2_train,upd.opt,scores_wapls,inf_train_wapls,result1,wapls_val[[5]])
        names(results)[[7]]<-"inferred train.set.val"
        }  
   if (val=="10-cross")
        {
        wapls_val<-tencross.wapls(train_set,train_env)
        result1[,5]<-wapls_val[[1]]
        result1[,6]<-wapls_val[[2]]
        result1[,7]<-wapls_val[[3]]
        result1[,8]<-wapls_val[[4]]
        results<-list(spec_n_train,n2_train,upd.opt,scores_wapls,inf_train_wapls,result1,wapls_val[[5]],wapls_val[[6]])
        names(results)[[7]]<-"mean(inferred train.set).val"
        names(results)[[8]]<-"sd(inferred train.set).val"
        }
   if (val=="boot")
        {
        wapls_val<-boot1.wapls(train_set,train_env)
        result1[,5]<-wapls_val[[1]]
        result1[,6]<-wapls_val[[2]]
        result1[,7]<-wapls_val[[3]]
        result1[,8]<-wapls_val[[4]]
        results<-list(spec_n_train,n2_train,upd.opt,scores_wapls,inf_train_wapls,result1,wapls_val[[5]],wapls_val[[6]],wapls_val[[7]],wapls_val[[8]])
        names(results)[[7]]<-"mean(inferred train.set).val"
        names(results)[[8]]<-"s1 (boot)"
        names(results)[[9]]<-"s2 (boot)"
        names(results)[[10]]<-"sd(inferred train.set).val"
        }
    
    names(results)[[1]]<-"species in train.set"
    names(results)[[2]]<-"N2 train.set"
    names(results)[[3]]<-"updated opt."
    names(results)[[4]]<-"sample scores"
    names(results)[[5]]<-"inferred train.set"
    names(results)[[6]]<-"performance"  
}
   
if (dat==3)
    {
    if (val=="none")
        {
        results<-list(spec_n_train,n2_train,upd.opt,scores_wapls,inf_train_wapls,result1,spec_n_test,n_test.train,n2_test,wapls_run2[[6]])
        names(results)[[7]]<-"species in core.samples"
        names(results)[[8]]<-"n species core.samples in train.set"
        names(results)[[9]]<-"N2 in core.samples"      
        names(results)[[10]]<-"reconstruction_core.samples"
        }
    
    if (val=="loo")
        {
        wapls_val<-loo2.wapls(train_set,train_env,test_set)
        result1[,5]<-wapls_val[[1]]
        result1[,6]<-wapls_val[[2]]
        result1[,7]<-wapls_val[[3]]
        result1[,8]<-wapls_val[[4]]
        results<-list(spec_n_train,n2_train,upd.opt,scores_wapls,inf_train_wapls,result1,wapls_val[[5]],spec_n_test,n_test.train,n2_test,wapls_run2[[6]],wapls_val[[6]],wapls_val[[7]])
        names(results)[[7]]<-"inferred train.set.val"
        names(results)[[8]]<-"species in core.samples"
        names(results)[[9]]<-"n species core.samples in train.set"
        names(results)[[10]]<-"N2 in core.samples"      
        names(results)[[11]]<-"reconstruction_core.samples"
        names(results)[[12]]<-"mean(reconstruction_core.samples).val"
        names(results)[[13]]<-"sd(reconstruction_core.samples).val"
        }
    
    if (val=="10-cross")
        {
        wapls_val<-tencross.wapls(train_set,train_env)
        result1[,5]<-wapls_val[[1]]
        result1[,6]<-wapls_val[[2]]
        result1[,7]<-wapls_val[[3]]
        result1[,8]<-wapls_val[[4]]
        results<-list(spec_n_train,n2_train,upd.opt,scores_wapls,inf_train_wapls,result1,wapls_val[[5]],wapls_val[[6]],spec_n_test,n_test.train,n2_test,wapls_run2[[6]])
        names(results)[[7]]<-"mean(inferred train.set).val"
        names(results)[[8]]<-"sd(inferred train.set).val"
        names(results)[[9]]<-"species in core.samples"
        names(results)[[10]]<-"n species core.samples in train.set"
        names(results)[[11]]<-"N2 in core.samples"      
        names(results)[[12]]<-"reconstruction_core.samples"
        }
     if (val=="boot")
        {
        wapls_val<-boot2.wapls(train_set,train_env,test_set)
        result1[,5]<-wapls_val[[12]]
        result1[,6]<-wapls_val[[8]]
        result1[,7]<-wapls_val[[9]]
        result1[,8]<-wapls_val[[11]]
        results<-list(spec_n_train,n2_train,upd.opt,scores_wapls,inf_train_wapls,result1,wapls_val[[5]],wapls_val[[6]],wapls_val[[7]],wapls_val[[13]],spec_n_test,n_test.train,n2_test,wapls_run2[[6]],wapls_val[[4]],wapls_val[[1]])
        names(results)[[7]]<-"mean(inferred train.set).val "
        names(results)[[8]]<-"s1 (boot)"
        names(results)[[9]]<-"s2 (boot)"
        names(results)[[10]]<-"sd(inferred train.set).val"
        names(results)[[11]]<-"species in core.samples"
        names(results)[[12]]<-"n species core.samples in train.set"
        names(results)[[13]]<-"N2 in core.samples"      
        names(results)[[14]]<-"reconstruction_core.samples"
        names(results)[[15]]<-"mean(reconstruction_core.samples).val"
        names(results)[[16]]<-"sd(reconstruction_core.samples).val"  
        }
    names(results)[[1]]<-"species in train.set"
    names(results)[[2]]<-"N2 train.set"
    names(results)[[3]]<-"updated opt."
    names(results)[[4]]<-"sample scores"
    names(results)[[5]]<-"inferred train.set"
    names(results)[[6]]<-"performance"  
    }

if (out=="TRUE")
{
if (val!="none")
    {
    cat("",fill=TRUE)
    if (val=="boot")
    {
    s.1.2<-matrix(ncol=2,nrow=n_comp,dimnames=list(nam1,c("s1","s2")))
    s.1.2[,1]<-wapls_val[[6]]
    s.1.2[,2]<-wapls_val[[7]]
    print(s.1.2)
    }
    cat("",fill=TRUE)
    print(round(result1[,5:8],7))
    cat("",fill=TRUE)
    cat("",fill=TRUE)
    cat("  RMSEP change [%]",fill=TRUE)
    cat("",fill=TRUE)
    for (i in 1:(comp-1))
        cat(paste("comp",i," : comp",i+1,":  ",round((result1[i,8]-result1[(i+1),8])*100/result1[i,8],2),sep=""),fill=TRUE)
    
    }
}   
########################################################################plots###############################################################################

if(d.plot==TRUE)
    {
     error1<-result1[order(result1[,4]),]
     comp.p<-error1[1,9]
     inf<-(inf_train_wapls[,comp.p])
     if (val=="none")
        {   
        par(mfrow=c(1,2))
        x1<-range(train_env)[1]-range(train_env)[1]*0.5
        x2<-range(train_env)[2]+range(train_env)[2]*0.2
        plot(inf~train_env,xlab="observed",ylab=paste("inferred (train_set)",nam1[comp.p]),col="gray",pch=19,xlim=c(x1,x2),ylim=c(x1,x2))
        abline(0,1,col="red",lty=2)
        panel.smooth(train_env,inf,col.smooth="blue",lwd=2)
        dif.plot<-(inf-train_env)
        plot(dif.plot~train_env,col="gray",pch=19,xlim=c(x1,x2),xlab="observed",ylab="residuals (inferred-measured)")
        abline(0,0,col="red",lty=2)
        panel.smooth(train_env,dif.plot,col.smooth="blue",lwd=2)
        }
    if (val!="none")
        {   
        error1<-result1[order(result1[,8]),]
        comp.p<-error1[1,9]
        inf<-inf_train_wapls[,comp.p]
        inf2<-wapls_val[[5]][,comp.p]
        par(mfrow=c(2,2))
        x1<-range(train_env)[1]-range(train_env)[1]*0.5
        x2<-range(train_env)[2]+range(train_env)[2]*0.2
        plot(inf~train_env,xlab="observed",ylab=paste("inferred (train_set)",nam1[comp.p]),col="gray",pch=19,xlim=c(x1,x2),ylim=c(x1,x2))
        abline(0,1,col="red",lty=2)
        panel.smooth(train_env,inf,col.smooth="blue",lwd=2)
        dif.plot<-(inf-train_env)
        plot(dif.plot~train_env,col="gray",pch=19,xlim=c(x1,x2),xlab="observed",ylab="residuals (inferred-measured)")
        abline(0,0,col="red",lty=2)
        panel.smooth(train_env,dif.plot,col.smooth="blue",lwd=2)
        plot(inf2~train_env,xlab="observed",ylab=paste("inferred (",val1,") (train_set)",nam1[comp.p]),col="gray",pch=19,xlim=c(x1,x2),ylim=c(x1,x2))
        abline(0,1,col="red",lty=2)
        panel.smooth(train_env,inf2,col.smooth="blue",lwd=2)
        dif.plot<-(inf2-train_env)
        plot(dif.plot~train_env,col="gray",pch=19,xlim=c(x1,x2),xlab="observed",ylab="residuals (cross.validation-measured)")
        abline(0,0,col="red",lty=2)
        panel.smooth(train_env,dif.plot,col.smooth="blue",lwd=2)
        }
    }












#############################################################################################################################################################     
results
}