
pom<-function(...,d.plot=TRUE,scale=TRUE,datatype=TRUE,min.occ=4,set.zero=FALSE,not.av=c("zero","lower","max"),val="loo",class=c(1,5,10,30,60),out=TRUE) 
{
    data<-list(...)
    train_set <- as.matrix(data[[1]])
    train_env <- as.matrix(data[[2]])
    #train_set <- as.matrix(train_set.MV)
    #train_env <- as.matrix(train_env.MV)
    env<-NA
    class.l<-length(class)
    if (missing(not.av)) 
        not.av<-"max"
        
    if (length(data) == 3) 
        test_set <- data[[3]]
    if (length(data) <= 2) 
        test_set <- NA
    train_set<-train_set[,colSums(train_set)!=0]
    n<-min.occ
    
    
  
    
    
    
    if (scale==TRUE)
    {
    if (length(data)==3)
        {test_set<-data[[3]]
        if (dim(test_set)[1]>1)
            test_set<-as.data.frame(t(apply(test_set,1,function(x)x*100/sum(x,na.rm=TRUE)))) 
        }

    train_set<-as.data.frame(t(apply(train_set,1,function(x)x*100/sum(x,na.rm=TRUE))))
    } 

      train_set.n<-train_set
    col<-dim(train_set.n)[2]
    row<-dim(train_set.n)[1] 
  
    nspec<-rep(NA,col)
    prob_spec<-as.list(nspec)
    
    env_mi<-min(train_env[,1])
    env_ma<-max(train_env[,1])
    test_tp<-NA
   
    delta<-(env_ma-env_mi)/3000
    for (i in 1:3000) test_tp[i]<-env_mi+i*delta
    train_set.pa<-train_set.n
    
    for (i in 1:col)
    {    
      fit.glm<-"NA"
      if (datatype==TRUE) 
        {
            spec.10<-train_set[,i]+10
            spec.pa<-spec.10
            for (cl in 1:class.l)
                spec.pa=ifelse(spec.pa>class[class.l-(cl-1)]+10,(class.l+1)-(cl-1),spec.pa)   
        
            spec.pa=ifelse(spec.pa>10,1,spec.pa)   
            spec.pa=ifelse(spec.pa>9,0,spec.pa)
            train_set.pa[,i]<-spec.pa
            names(prob_spec)[[i]]<-names(train_set.n)[i] 
        }
      for (ln in (length(as.numeric(levels(as.factor(spec.pa))))):2)
        {
          if (length(spec.pa[spec.pa==as.numeric(levels(as.factor(spec.pa)))[ln]])<=4)
            {
                if (set.zero==TRUE)
                    spec.pa[spec.pa==as.numeric(levels(as.factor(spec.pa)))[ln]]<-0
                if (set.zero==FALSE)
                    spec.pa[spec.pa==as.numeric(levels(as.factor(spec.pa)))[ln]]<-as.numeric(levels(as.factor(spec.pa)))[ln-1]
            }
        }
      if(length(levels(as.factor(spec.pa)))>2)
        {
            fit2.polr<-"NA"
            env<-train_env[,1]
            res<-matrix(ncol=length(levels(as.factor(spec.pa))),nrow=3000)  
            spec.pa0<-ifelse(spec.pa>1,1,spec.pa) 
            fit.glm<-glm(spec.pa0~env+I(env^2),quasibinomial)
            env<-test_tp
            fit.pre<-predict(fit.glm,newdata=as.data.frame(env),type="response")
            lev<-levels(as.factor(spec.pa))
            colnames(res)<-lev
            res[,1]<-(1-fit.pre) 
            res<-as.data.frame(res) 
            for (lev.run in 2:length(levels(as.factor(spec.pa))))
                {
                    spec.pa0<-ifelse(spec.pa==as.numeric(levels(as.factor(spec.pa))[lev.run]),1,0) 
                    env<-train_env[,1]
                    fit.glm<-glm(spec.pa0~env+I(env^2),quasibinomial)
                    env<-test_tp
                    fit.pre<-predict(fit.glm,newdata=as.data.frame(env),type="response")
                    res[,lev.run]<-(fit.pre) 
                }
            prob_spec[[i]]<-res 
        }
    if(length(levels(as.factor(spec.pa)))==2)
        {
            env<-train_env[,1]
            lev<-levels(as.factor(spec.pa))
            spec.pa<-ifelse(spec.pa>1,1,spec.pa)
            fit.glm<-glm(spec.pa~env+I(env^2),quasibinomial)
            
            env<-test_tp
            fit.pre<-predict(fit.glm,newdata=as.data.frame(env),type="response")
            res<-matrix(ncol=2,nrow=3000)
            colnames(res)<-lev
            res[,1]<-(1-fit.pre)
            res[,2]<-fit.pre
            res<-as.data.frame(res)
            prob_spec[[i]]<-res
        }

    }

############################ predicton on the train_set ###############################################
    end_tp<-NA
    for (j in 1:row)
        {
            end_tp_try<-rep(1,3000)
            for (i in 1:col) 
                {
                    if(is.na(sum((prob_spec[[i]])))==FALSE)
                        {
                         {
                            k<-paste(train_set.pa[j,i])
           #if (dim(prob_spec[[i]])[2]>=(train_set.pa[j,i]+1))
                            if(k%in%names(prob_spec[[i]]))
                                 end_tp_try<-end_tp_try*prob_spec[[i]][,k]
                          
                          
                          if (not.av=="max")
                            {
                                if (as.numeric(k)>max(as.numeric(names(prob_spec[[i]]))))
                                    { 
                                        if (dim(prob_spec[[i]])[2]>2)
                                            end_tp_try<-end_tp_try*rowSums(prob_spec[[i]][,2:dim(prob_spec[[i]])[2]])
                                        if (dim(prob_spec[[i]])[2]==2)
                                            end_tp_try<-end_tp_try*prob_spec[[i]][,dim(prob_spec[[i]])[2]]
                                    }
                            }
                          if (not.av=="lower")  
                                if (as.numeric(k)>max(as.numeric(names(prob_spec[[i]]))))
                                        end_tp_try<-end_tp_try*prob_spec[[i]][,dim(prob_spec[[i]])[2]]
                          if (not.av=="zero") 
                                if (as.numeric(k)>max(as.numeric(names(prob_spec[[i]]))))
                                        end_tp_try<-end_tp_try
                        }
                       }
                }
            tpdat<-cbind(end_tp_try,test_tp)
            tpdat<-tpdat[order(tpdat[,1]),]
            end_tp[j]<-tpdat[3000,2]
        }
    RMSE<-sqrt((sum((train_env-end_tp)^2))/length(train_env))
    R2<-summary(lm(train_env~end_tp))$r.squared
    max_error<-max(train_env-end_tp)
    mean_error<-mean(train_env-end_tp)
if (out==TRUE)
   {
    cat("",fill=TRUE)
    cat("",fill=TRUE)
    cat("                   transfer function",fill=TRUE)
    cat("",fill=TRUE)
    cat("",fill=TRUE)
    cat("type        = pom",fill=TRUE)
    cat("n samples   =",dim(train_set)[1],fill=TRUE)
    cat("n species   =",dim(train_set)[2],fill=TRUE)
    cat("val.-method = leave-one-out",fill=TRUE) 
    cat("",fill=TRUE)
    cat("R2    = ",round(R2,4),"          Mean-error = ",round(mean_error,4),fill=TRUE)
    cat("RMSE  = ",round(RMSE,4),"          Max-error  = ",round(max_error,4),fill=TRUE)
    cat("",fill=TRUE)
    if (val=="loo")
        { 
            cat("                     * please wait *",fill=TRUE) 
            cat("",fill=TRUE)
            cat("",fill=TRUE)
        }
    flush.console()
   }
###################################################          test_set    ################################################
    env<-NA
    #prob_spec.t<-NA
    if (length(data) == 3) 
    {
        end_tp_test<-NA
        env<-train_env[,1]
        k<-colnames(test_set) %in% colnames(train_set)       
        k<-colnames(test_set)[k]
        test_set.n<-test_set[,k]
        train_set.n<-train_set[k]
        col<-dim(train_set.n)[2]
        row<-dim(train_set.n)[1] 
        nspec<-rep(NA,col)
        prob_spec.t<-as.list(nspec)
        test_set.pa<-test_set.n
        train_set1.pa<-train_set.n
        for (i in 1:col)
            {    
            fit.glm<-"NA"
            if (datatype==TRUE) 
                 {
                    spec.10<-train_set.n[,i]+10
                    
                    spec.pa<-spec.10
                    for (cl in 1:class.l)
                        spec.pa=ifelse(spec.pa>class[class.l-(cl-1)]+10,(class.l+1)-(cl-1),spec.pa)   
                    spec.pa=ifelse(spec.pa>10,1,spec.pa)   
                    spec.pa=ifelse(spec.pa>9,0,spec.pa)
                    train_set1.pa[,i]<-spec.pa
                    names(prob_spec.t)[[i]]<-names(train_set.n)[i] 
                 }
            for (ln in (length(as.numeric(levels(as.factor(spec.pa))))):2)
                {
                    if (length(spec.pa[spec.pa==as.numeric(levels(as.factor(spec.pa)))[ln]])<=4)
                        {
                            if (set.zero==TRUE)
                               spec.pa[spec.pa==as.numeric(levels(as.factor(spec.pa)))[ln]]<-0
                            if (set.zero==FALSE)
                                spec.pa[spec.pa==as.numeric(levels(as.factor(spec.pa)))[ln]]<-as.numeric(levels(as.factor(spec.pa)))[ln-1]
                        }
                }
            if(length(levels(as.factor(spec.pa)))>2)
                {
                    fit2.polr<-"NA"
                    env<-NA
                    env<-train_env[,1]
                    res<-matrix(ncol=length(levels(as.factor(spec.pa))),nrow=3000)  
                    spec.pa0<-ifelse(spec.pa>1,1,spec.pa) 
                    fit.glm<-glm(spec.pa0~env+I(env^2),quasibinomial)
                    env<-test_tp
                    fit.pre<-predict(fit.glm,newdata=as.data.frame(env),type="response")
                    lev<-levels(as.factor(spec.pa))
                    colnames(res)<-lev
                    res[,1]<-(1-fit.pre) 
                    res<-as.data.frame(res) 
                    for (lev.run in 2:length(levels(as.factor(spec.pa))))
                        {
                            spec.pa0<-ifelse(spec.pa==as.numeric(levels(as.factor(spec.pa))[lev.run]),1,0) 
                            env<-NA
                            env<-train_env[,1]
                            fit.glm<-glm(spec.pa0~env+I(env^2),quasibinomial)
                            env<-test_tp
                            fit.pre<-predict(fit.glm,newdata=as.data.frame(env),type="response")
                            res[,lev.run]<-(fit.pre) 
                        }
                    prob_spec.t[[i]]<-res 
                }
            if(length(levels(as.factor(spec.pa)))==2)
                {
                    env<-train_env[,1]
                    lev<-levels(as.factor(spec.pa))
                    spec.pa<-ifelse(spec.pa>1,1,spec.pa)
                    fit.glm<-glm(spec.pa~env+I(env^2),quasibinomial)
                    env<-test_tp
                    fit.pre<-predict(fit.glm,newdata=as.data.frame(env),type="response")
                    res<-matrix(ncol=2,nrow=3000)
                    colnames(res)<-lev
                    res[,1]<-(1-fit.pre)
                    res[,2]<-fit.pre
                    res<-as.data.frame(res)
                    prob_spec.t[[i]]<-res
             }
    }

#######  prediction for test_set   ######
 
  
  
  if (datatype==TRUE) 
                 {
                  for (i in 1:dim(test_set.n)[2])
                  
                    {spec.10<-test_set.n[,i]+10
                    
                    spec.pa<-spec.10
                    for (cl in 1:class.l)
                         spec.pa=ifelse(spec.pa>class[class.l-(cl-1)]+10,(class.l+1)-(cl-1),spec.pa)   
                    spec.pa=ifelse(spec.pa>10,1,spec.pa)   
                    spec.pa=ifelse(spec.pa>9,0,spec.pa)
         
                    test_set.pa[,i]<-spec.pa}
                    #names(prob_spec.t)[[i]]<-names(train_set.n)[i] 
                 }
  
  
  
  
  for (j in 1:dim(test_set.n)[1])
        {
            end_tp_try<-rep(1,3000)
            for (i in 1:dim(test_set.n)[2]) 
                {
                   if(is.na(sum((prob_spec.t[[i]])))==FALSE)
                     {
                       k<-paste(test_set.pa[j,i])
                        #if (dim(prob_spec[[i]])[2]>=(train_set.pa[j,i]+1))
                       if(k%in%names(prob_spec.t[[i]]))
                            end_tp_try<-end_tp_try*prob_spec.t[[i]][,k]
                       if (not.av=="max")
                            {
                                if (as.numeric(k)>max(as.numeric(names(prob_spec.t[[i]])),na.rm = TRUE))
                                    { 
                                        if (dim(prob_spec.t[[i]])[2]>2)
                                            end_tp_try<-end_tp_try*rowSums(prob_spec.t[[i]][,2:dim(prob_spec.t[[i]])[2]])
                                        if (dim(prob_spec.t[[i]])[2]==2)
                                            end_tp_try<-end_tp_try*prob_spec.t[[i]][,dim(prob_spec.t[[i]])[2]]
                                    }
                            }
                          if (not.av=="lower")  
                                if (as.numeric(k)>max(as.numeric(names(prob_spec.t[[i]]))))
                                        end_tp_try<-end_tp_try*prob_spec.t[[i]][,dim(prob_spec.t[[i]])[2]]
                          if (not.av=="zero") 
                                if (as.numeric(k)>max(as.numeric(names(prob_spec.t[[i]]))))
                                        end_tp_try<-end_tp_try
                    }
                }
            #plot(end_tp_try[,1]~test_tp)
            tpdat<-cbind(end_tp_try,test_tp)
            tpdat<-tpdat[order(tpdat[,1]),]
            end_tp_test[j]<-tpdat[3000,2]
        }

}

#######################################################    LOO           ################################################

if (datatype==TRUE) 
        {
            
            spec.10<-train_set[,i]+10
            spec.pa<-spec.10
            for (cl in 1:class.l)
                spec.pa=ifelse(spec.pa>class[class.l-(cl-1)]+10,(class.l+1)-(cl-1),spec.pa)   
            spec.pa=ifelse(spec.pa>10,1,spec.pa)   
            spec.pa=ifelse(spec.pa>9,0,spec.pa)
            train_set.pa[,i]<-spec.pa
            names(prob_spec)[[i]]<-names(train_set.n)[i] 
        }

if (val=="loo")
 {
    
    
    
    end_tp1<-NA
    env<-NA
    for (loo in 1:dim(train_set)[1])
        {
        train_set.paloo<-train_set.pa[-loo,]
        train_env.loo<-train_env[-loo,]
        col<-dim(train_set.paloo)[2]
        nspec<-rep(NA,col)
        prob_spec.l<-as.list(nspec)
        
        for (i in 1:col)
            {    
                spec.pa<-train_set.paloo[,i]
                for (ln in (length(as.numeric(levels(as.factor(spec.pa))))):2)
                  {
                    if (length(spec.pa[spec.pa==as.numeric(levels(as.factor(spec.pa)))[ln]])<=4)
                        {
                            if (set.zero==TRUE)
                                spec.pa[spec.pa==as.numeric(levels(as.factor(spec.pa)))[ln]]<-0
                            if (set.zero==FALSE)
                                 spec.pa[spec.pa==as.numeric(levels(as.factor(spec.pa)))[ln]]<-as.numeric(levels(as.factor(spec.pa)))[ln-1]
                        }
                    }
        if(length(levels(as.factor(spec.pa)))>2)
            {
                fit2.polr<-"NA"
                env<-train_env.loo#[,1]
                res<-matrix(ncol=length(levels(as.factor(spec.pa))),nrow=3000)  
                spec.pa0<-ifelse(spec.pa>1,1,spec.pa) 
                fit.glm<-glm(spec.pa0~env+I(env^2),quasibinomial)
                env<-test_tp
                fit.pre<-predict(fit.glm,newdata=as.data.frame(env),type="response")
                lev<-levels(as.factor(spec.pa))
                colnames(res)<-lev
                res[,1]<-(1-fit.pre) 
                res<-as.data.frame(res) 
                for (lev.run in 2:length(levels(as.factor(spec.pa))))
                    {
                    spec.pa0<-ifelse(spec.pa==as.numeric(levels(as.factor(spec.pa))[lev.run]),1,0) 
                    env<-train_env.loo
                    fit.glm<-glm(spec.pa0~env+I(env^2),quasibinomial)                    
                    env<-test_tp
                    fit.pre<-predict(fit.glm,newdata=as.data.frame(env),type="response")
                    res[,lev.run]<-(fit.pre) 
                    }        
                prob_spec.l[[i]]<-res       
            }
        if(length(levels(as.factor(spec.pa)))==2)
            {
                env<-train_env.loo
                lev<-levels(as.factor(spec.pa))
                spec.pa<-ifelse(spec.pa>1,1,spec.pa)
                fit.glm<-glm(spec.pa~env+I(env^2),quasibinomial)               
                env<-test_tp
                fit.pre<-predict(fit.glm,newdata=as.data.frame(env),type="response")
                res<-matrix(ncol=2,nrow=3000)
                colnames(res)<-lev
                res[,1]<-(1-fit.pre)
                res[,2]<-fit.pre
                res<-as.data.frame(res)
                prob_spec.l[[i]]<-res
            }
    }

############################ predicton on the train_set ###############################################


    end_tp_try<-rep(1,3000)
    for (i in 1:col) 
        {
            if(is.na(sum((prob_spec.l[[i]])))==FALSE)
                {
                     k<-paste(train_set.pa[loo,i])
                    # if (dim(prob_spec.l[[i]])[2]>=(train_set.pa[loo,i]+1))
                        if(k%in%names(prob_spec.l[[i]]))
                                end_tp_try<-end_tp_try*prob_spec.l[[i]][,k]
                     if (not.av=="max")
                            {
                                if (as.numeric(k)>max(as.numeric(names(prob_spec.l[[i]]))))
                                    { 
                                        if (dim(prob_spec.l[[i]])[2]>2)
                                            end_tp_try<-end_tp_try*rowSums(prob_spec.l[[i]][,2:dim(prob_spec.l[[i]])[2]])
                                        if (dim(prob_spec[[i]])[2]==2)
                                            end_tp_try<-end_tp_try*prob_spec.l[[i]][,dim(prob_spec.l[[i]])[2]]
                                    }
                            }
                     if (not.av=="lower")  
                            if (as.numeric(k)>max(as.numeric(names(prob_spec.l[[i]]))))
                                end_tp_try<-end_tp_try*prob_spec.l[[i]][,dim(prob_spec.l[[i]])[2]]
                     if (not.av=="zero") 
                            if (as.numeric(k)>max(as.numeric(names(prob_spec.l[[i]]))))
                                end_tp_try<-end_tp_try
                }
        

       
        tpdat<-cbind(end_tp_try,test_tp)
        tpdat<-tpdat[order(tpdat[,1]),]
        end_tp1[loo]<-tpdat[3000,2]
        }

      }
   
}
###############                                  loo for test_set                                                                                        
                                  
if (length(data)==3 & val=="loo")                   
    {
       end_tp_test.m<-matrix(ncol=dim(train_set)[1],nrow=dim(test_set)[1])
       for (lo.t in 1:dim(train_set)[1]) 
        {
        train_set.l<-train_set[-lo.t,]
        train_env.l<-as.matrix(train_env[-lo.t,])
        train_set.l<-train_set.l[,colSums(train_set.l)!=0]
        
        env<-train_env.l[,1]
        k<-colnames(test_set) %in% colnames(train_set.l)       
        k<-colnames(test_set)[k]
        test_set.n.loo<-test_set[,k]
        train_set.n.loo<-train_set.l[k]
        col<-dim(train_set.n.loo)[2]
        row<-dim(train_set.n.loo)[1] 
        nspec<-rep(NA,col)
        prob_spec.t<-as.list(nspec)
        test_set.pa<-test_set.n.loo
        train_set1.pa.loo<-train_set.n.loo
        for (i in 1:col)
            {    
            fit.glm<-"NA"
            if (datatype==TRUE) 
                 {
                    spec.10<-train_set.n.loo[,i]+10
                    spec.pa<-spec.10
                    for (cl in 1:class.l)
                        spec.pa=ifelse(spec.pa>class[class.l-(cl-1)]+10,(class.l+1)-(cl-1),spec.pa)   
                    spec.pa=ifelse(spec.pa>10,1,spec.pa)   
                    spec.pa=ifelse(spec.pa>9,0,spec.pa)
                    train_set1.pa.loo[,i]<-spec.pa
                    names(prob_spec.t)[[i]]<-names(train_set.n.loo)[i] 
                 }
            for (ln in (length(as.numeric(levels(as.factor(spec.pa))))):2)
                {
                    if (length(spec.pa[spec.pa==as.numeric(levels(as.factor(spec.pa)))[ln]])<=4)
                        {
                            if (set.zero==TRUE)
                               spec.pa[spec.pa==as.numeric(levels(as.factor(spec.pa)))[ln]]<-0
                            if (set.zero==FALSE)
                                spec.pa[spec.pa==as.numeric(levels(as.factor(spec.pa)))[ln]]<-as.numeric(levels(as.factor(spec.pa)))[ln-1]
                        }
                }
            if(length(levels(as.factor(spec.pa)))>2)
                {
                    fit2.polr<-"NA"
                    env<-NA
                    env<-train_env.l[,1]
                    res<-matrix(ncol=length(levels(as.factor(spec.pa))),nrow=3000)  
                    spec.pa0<-ifelse(spec.pa>1,1,spec.pa) 
                    fit.glm<-glm(spec.pa0~env+I(env^2),quasibinomial)
                    env<-test_tp
                    fit.pre<-predict(fit.glm,newdata=as.data.frame(env),type="response")
                    lev<-levels(as.factor(spec.pa))
                    colnames(res)<-lev
                    res[,1]<-(1-fit.pre) 
                    res<-as.data.frame(res) 
                    for (lev.run in 2:length(levels(as.factor(spec.pa))))
                        {
                            spec.pa0<-ifelse(spec.pa==as.numeric(levels(as.factor(spec.pa))[lev.run]),1,0) 
                            env<-NA
                            env<-train_env.l[,1]
                            fit.glm<-glm(spec.pa0~env+I(env^2),quasibinomial)
                            env<-test_tp
                            fit.pre<-predict(fit.glm,newdata=as.data.frame(env),type="response")
                            res[,lev.run]<-(fit.pre) 
                        }
                    prob_spec.t[[i]]<-res 
                }
            if(length(levels(as.factor(spec.pa)))==2)
                {
                    env<-train_env.l[,1]
                    lev<-levels(as.factor(spec.pa))
                    spec.pa<-ifelse(spec.pa>1,1,spec.pa)
                    fit.glm<-glm(spec.pa~env+I(env^2),quasibinomial)
                    env<-test_tp
                    fit.pre<-predict(fit.glm,newdata=as.data.frame(env),type="response")
                    res<-matrix(ncol=2,nrow=3000)
                    colnames(res)<-lev
                    res[,1]<-(1-fit.pre)
                    res[,2]<-fit.pre
                    res<-as.data.frame(res)
                    prob_spec.t[[i]]<-res
             }
    }

#######  prediction for test_set   ######
 
  
  
  if (datatype==TRUE) 
                 {
                  for (i in 1:dim(test_set.n.loo)[2])
                  
                    {spec.10<-test_set.n.loo[,i]+10
                        spec.pa<-spec.10
                        for (cl in 1:class.l)
                            spec.pa=ifelse(spec.pa>class[class.l-(cl-1)]+10,(class.l+1)-(cl-1),spec.pa)   
                        spec.pa=ifelse(spec.pa>10,1,spec.pa)   
                        spec.pa=ifelse(spec.pa>9,0,spec.pa)
                    test_set.pa[,i]<-spec.pa}
                    #names(prob_spec.t)[[i]]<-names(train_set.n.loo)[i] 
                 }
  
  
  
  
  for (j in 1:dim(test_set.n.loo)[1])
        {
            end_tp_try<-rep(1,3000)
            for (i in 1:dim(test_set.n.loo)[2]) 
                {
                   if(is.na(sum((prob_spec.t[[i]])))==FALSE)
                     {
                       k<-paste(test_set.pa[j,i])
                        #if (dim(prob_spec[[i]])[2]>=(train_set.pa[j,i]+1))
                       if(k%in%names(prob_spec.t[[i]]))
                            end_tp_try<-end_tp_try*prob_spec.t[[i]][,k]
                       if (not.av=="max")
                            {
                                if (as.numeric(k)>max(as.numeric(names(prob_spec.t[[i]])),na.rm = TRUE))
                                    { 
                                        if (dim(prob_spec.t[[i]])[2]>2)
                                            end_tp_try<-end_tp_try*rowSums(prob_spec.t[[i]][,2:dim(prob_spec.t[[i]])[2]])
                                        if (dim(prob_spec.t[[i]])[2]==2)
                                            end_tp_try<-end_tp_try*prob_spec.t[[i]][,dim(prob_spec.t[[i]])[2]]
                                    }
                            }
                          if (not.av=="lower")  
                                if (as.numeric(k)>max(as.numeric(names(prob_spec.t[[i]]))))
                                        end_tp_try<-end_tp_try*prob_spec.t[[i]][,dim(prob_spec.t[[i]])[2]]
                          if (not.av=="zero") 
                                if (as.numeric(k)>max(as.numeric(names(prob_spec.t[[i]]))))
                                        end_tp_try<-end_tp_try
                    }
                }
            #plot(end_tp_try[,1]~test_tp)
            tpdat<-cbind(end_tp_try,test_tp)
            tpdat<-tpdat[order(tpdat[,1]),]
            end_tp_test.m[j,lo.t]<-tpdat[3000,2]
        }
     end_tp_test.mean<-NA
     end_tp_test.sd<-NA
     for (k in 1:dim(test_set)[1])
        {
            end_tp_test.mean[k]<-mean(end_tp_test.m[k,])
            end_tp_test.sd[k]<-sd(end_tp_test.m[k,])
        }
    }
}   

  
  
  
  
  
  
  

##################################         end loo                                                                                

train_env<-train_env[,1]
if (val=="loo")
    {
        error_m<-matrix(nrow=1,ncol=8,dimnames=list("pom",c("RMSE","R2","Ave_Bias","Max_Bias","X_R2","X_Ave_Bias","X_Max_Bias","RMSEP")))
        RMSEP<-sqrt((sum((train_env-end_tp1)^2))/length(train_env))
        R2_cv<-summary(lm(train_env~end_tp1))$r.squared
        max_error_cv<-max(train_env-end_tp1)
        mean_error_cv<-mean(train_env-end_tp1)
      
 if (out==TRUE)
   {
        cat("R2.c     = ",R2_cv,"          Mean-error.c     = ",round(mean_error_cv,4),fill=TRUE)
        cat("RMSEP    = ",RMSEP,"          Max-error.c      = ",round(max_error_cv,4),fill=TRUE)
        cat("",fill=TRUE)
        cat("",fill=TRUE)
      
   }   
        error_m[1,1]<-RMSE
        error_m[1,2]<-R2
        error_m[1,3]<-mean_error
        error_m[1,4]<-max_error
        error_m[1,5]<-R2_cv
        error_m[1,6]<-mean_error_cv
        error_m[1,7]<-max_error_cv
        error_m[1,8]<-RMSEP
        
        if (length(data)==2)
           {
           results<-list(end_tp,error_m,prob_spec,end_tp1)          
           names(results)[[1]]<-"inferred train.set"
           names(results)[[2]]<-"performance"
           names(results)[[3]]<-"spec distribution"
           names(results)[[4]]<-"inferred train.set (loo)"
          
          
           }
        
        if (length(data)==3)
           {
           results<-list(end_tp,error_m,prob_spec,end_tp1,end_tp_test,end_tp_test.mean,end_tp_test.sd)
           names(results)[[1]]<-"inferred train.set"
           names(results)[[2]]<-"performance"
           names(results)[[3]]<-"spec distribution"
           names(results)[[4]]<-"inferred train.set (loo)"
           names(results)[[5]]<-"reconstruction_core.samples"
           names(results)[[6]]<-"mean (reconstruction_core.samples) (loo)"
           names(results)[[7]]<-"sd (reconstruction_core.samples) (loo)"
           }
}

if (val!="loo")
   {
    
    error_m<-matrix(nrow=1,ncol=8,dimnames=list("pom",c("RMSE","R2","Ave_Bias","Max_Bias","X_R2","X_Ave_Bias","X_Max_Bias","RMSEP")))
    error_m[1,1]<-RMSE
    error_m[1,2]<-R2
    error_m[1,3]<-mean_error
    error_m[1,4]<-max_error
    
    if (length(data)==2)
        {
        
        results<-list(end_tp,error_m,prob_spec)
        names(results)[[1]]<-"inferred train.set"
        names(results)[[2]]<-"performance"
        names(results)[[3]]<-"spec distribution"
        
        }
    
    if (length(data)==3)
        {
        results<-list(end_tp,error_m,prob_spec,end_tp_test)
        
        names(results)[[1]]<-"inferred train.set"
        names(results)[[2]]<-"performance"
        names(results)[[3]]<-"spec distribution"
        names(results)[[4]]<-"reconstruction_core.samples"
        
   }
   }



if (d.plot==TRUE)
    {
        if (val!="loo")
            {
                par(mfrow=c(1,2))
                x1<-range(train_env)[1]-range(train_env)[1]*0.5
                x2<-range(train_env)[2]+range(train_env)[2]*0.2
                plot(end_tp~train_env,xlab="observed",ylab="inferred (train_set)",col="gray",pch=19,xlim=c(x1,x2),ylim=c(x1,x2))
                abline(0,1,col="red",lty=2)
                panel.smooth(train_env,end_tp,col.smooth="blue",lwd=2)
                dif.plot<-(end_tp-train_env)
                plot(dif.plot~train_env,col="gray",pch=19,xlim=c(x1,x2),xlab="observed",ylab="residuals (inferred-measured)")
                abline(0,0,col="red",lty=2)
                panel.smooth(train_env,dif.plot,col.smooth="blue",lwd=2)
            }
       if (val=="loo")
            {
                par(mfrow=c(2,2))
                x1<-range(train_env)[1]-range(train_env)[1]*0.5
                x2<-range(train_env)[2]+range(train_env)[2]*0.2
                plot(end_tp~train_env,xlab="observed",ylab="inferred (train_set)",col="gray",pch=19,xlim=c(x1,x2),ylim=c(x1,x2))
                abline(0,1,col="red",lty=2)
                panel.smooth(train_env,end_tp,col.smooth="blue",lwd=2)
                dif.plot<-(end_tp-train_env)
                plot(dif.plot~train_env,col="gray",pch=19,xlim=c(x1,x2),xlab="observed",ylab="residuals (inferred-measured)")
                abline(0,0,col="red",lty=2)
                panel.smooth(train_env,dif.plot,col.smooth="blue",lwd=2)
            
                x1<-range(train_env)[1]-range(train_env)[1]*0.5
                x2<-range(train_env)[2]+range(train_env)[2]*0.2
                plot(end_tp1~train_env,xlab="observed",ylab="inferred (train_set)(loo)",col="gray",pch=19,xlim=c(x1,x2),ylim=c(x1,x2))
                abline(0,1,col="red",lty=2)
                panel.smooth(train_env,end_tp1,col.smooth="blue",lwd=2)
                dif.plot<-(end_tp1-train_env)
                plot(dif.plot~train_env,col="gray",pch=19,xlim=c(x1,x2),xlab="observed",ylab="residuals (inferred-measured)")
                abline(0,0,col="red",lty=2)
                panel.smooth(train_env,dif.plot,col.smooth="blue",lwd=2)
            
            
            }
    
    }
invisible(results)

}