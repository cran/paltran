

mw<-function(train_set,train_env,core_data,method=c("wapls","wa","pom"),comp=4,val=c("boot","loo","10-cross"),run=100,mwsize=c(20,40,60),dim=c(2,3,4),mw.type=c("dca","ca","cca","sample"),dist.m="euclidean",rmsep.incl=TRUE,env.trans=FALSE,spec.trans=FALSE,rplot=TRUE,drop.non.sig=FALSE,min.occ=1,scale=FALSE,dw=FALSE)

{
test_set<-core_data
train_set<-train_set[,colSums(train_set)!=0]
test_set<-test_set[,colSums(test_set)!=0]
ncomp<-comp



if (missing(dim)) 
        dim=2
if (missing(val)) 
        val <- "boot"
val<- match.arg(val)

if (missing(mw.type)) 
    mw.type<-"dca"

if (mw.type=="sample") 
    rplot<-FALSE
if (mw.type=="cca") 
    rplot<-FALSE




if (missing(method)) 
    method<-"wapls"  
    
if (method=="wa")
    ncomp=1

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



####################################################################################################################################################

best.wapls<-function(error,nro=4,rmsep.incl=TRUE)
    {
        note<-rep(0,nro)           
        if (nro<=2)
            {
            note[1]<-2
            note[2]<-1
            }             
        if (nro>=3)
            {
            note[1]<-2
            note[2]<-1
            note[3]<-0
            }
        if (nro>=2)    
                {
                error<-as.matrix(error)
                error[,10]<-rep(0,nro)
                #error[,9]<-c(1:nro)
                c<-abs(error[,6])
                error<-error[order(c),]
                error[,10]<-note
                    c<-error[,7]
                    error<-error[order(c),]
                    error[,10]<-error[,10]+note
                        c<-error[,5]
                        error<-error[order(c,decreasing = TRUE),]
                        error[,10]<-error[,10]+note
                if(rmsep.incl==TRUE)
                    {
                        c<-error[,8]        
                        error<-error[order(c,decreasing = FALSE),]
                        error[,10]<-error[,10]+note
                    }
                c<-error[,10]
                k<-error[,6]
                error<-error[order(c,k,decreasing=FALSE),]
                }  
    
    if (nro<=1)
       error[,9]<-c(1:nro)

    error

    }


#####################################################################################################################################################
row<-dim(train_set)[1]
col1<-dim(train_set)[2]
row1<-dim(test_set)[1]
number.mw_ds<-length(mwsize)
mw_env<-c(1,number.mw_ds)
mw_env.val<-c(1,number.mw_ds)
mw_env.sd<-c(1,number.mw_ds)
mw_env.end.val<-rep(NA,row1)
mw_env.end.sd<-rep(NA,row1)
mw_env.end<-rep(0,row1)
wapls_comp.mw<-matrix(nrow=(number.mw_ds),ncol=10)
nam1<-row.names(test_set)
mw_envstat.end<-matrix(NA,nrow=row1,ncol=10,dimnames=list(nam1,c("Rs","Ave_Bias","Max_Bias","RMSE","Rs.X","Ave_Bias.X","Max_Bias.X","RMSEP","comp","w.size")))

new<-merge(train_set,test_set,all=TRUE,sort=FALSE)                      # creating test and training_set with the same col
new<-as.data.frame(new)
row.names(new)<-c(row.names(train_set),row.names(test_set))
new<-new[,(1:col1)]
train_set<-new[(1:row),]
train_set<-apply(train_set,2,function(x){ifelse(is.na(x),0,x)})
test_set<-new[((row+1):dim(new)[1]),]
test_set<-apply(test_set,2,function(x){ifelse(is.na(x),0,x)})
new<-apply(new,2,function(x){ifelse(is.na(x),0,x)})
wapls_comp<-matrix(0,ncol=10,nrow=ncomp)


######################################################################################

library(vegan)
if (mw.type=="dca")
    {
    if (dw==TRUE)
        train.dca<-decorana(downweight(train_set)) 
    if (dw!=TRUE)
        train.dca<-decorana(train_set) 
    train_scores<-scores(train.dca)[,c(1:dim)]
    test.dca<-predict(train.dca,newdata=test_set,type="sites")  
    test_scores<-scores(test.dca)[,c(1:dim)]
    test_train.co<-rbind(train_scores,test_scores)
    dist<-vegdist(test_train.co,method=dist.m )
    dist<-as.matrix(dist)
    dist<-as.data.frame(dist)
    
    }

if (mw.type=="ca")
    {
    if (dw==TRUE)
        train.cca<-cca(downweight(train_set)) 
    if (dw!=TRUE)
        train.cca<-cca(train_set) 
    train_scores<-scores(train.cca,choices=c(1:dim),display="wa")
    test.cca<-predict(train.cca,newdata=test_set,type="wa") 
    test_scores<-scores(test.cca,choices=c(1:dim),display="sites")
    test_train.co<-rbind(train_scores,test_scores)
    dist<-vegdist(test_train.co,method=dist.m )
    dist<-as.matrix(dist)
    dist<-as.data.frame(dist)
    
    }
if (mw.type=="cca")
    {
    
        train.cca<-cca(downweight(train_set)~train_env[,1]) 
    if (dw!=TRUE)
        train.cca<-cca(train_set~train_env[,1])    
    train_scores<-scores(train.cca,choices=c(1:dim))$sites
    test.cca<-predict(train.cca,newdata=test_set,type="wa") 
    test_scores<-scores(test.cca)
    test_train.co<-rbind(as.matrix(train_scores[,1]),test_scores)
    dist<-vegdist(test_train.co,method=dist.m )
    dist<-as.matrix(dist)
    dist<-as.data.frame(dist)
    
    }


if (mw.type=="sample")
    {
    if (dw==TRUE)
        dist<-vegdist(downweight(new),method=dist.m)
    if (dw!=TRUE)
    dist<-vegdist(new,method=dist.m)
    dist<-as.matrix(dist)
    dist<-as.data.frame(dist)
   
    }

row<-dim(train_set)[1]
dist1<-dist[1:row,]

cat("",fill=TRUE)
cat("              moving window - transfer function",fill=TRUE)
cat("",fill=TRUE)
cat("",fill=TRUE)
cat("nearest neigbour method =",mw.type,fill=TRUE)
if (mw.type!="sample")
    cat("     (used dimension    =",dim,")",fill=TRUE)
cat("   distance measurement =",dist.m,fill=TRUE)
cat("   transfer function    =",method,fill=TRUE)
if (method=="wapls")
    cat("   number of components =",comp,fill=TRUE)
cat("   window size          =",mwsize,fill=TRUE)
cat("   val.-method          =",val,fill=TRUE)
if(val%in%c("boot","loo")) 
    cat("                 run =",run,fill=TRUE)
cat("",fill=TRUE)

cat("",fill=TRUE)
for(i in 1:row1)
    {
        for (k in 1:(number.mw_ds))
            {
                if (rplot==TRUE && mw.type!="cca")
                {
                    par(fig=c(0,1,0,0.9))
                    par(mar=c(5,5,0,2))
                    plot(train_scores[,1],train_scores[,2])
                    points(test_scores[,1],test_scores[,2],col="green",pch=19)
                    points(test_scores[c(1:i),1],test_scores[c(1:i),2],col="orange",pch=19) 
                }
                
                 
                
                
                
                
                nn.number<-mwsize[k]
                start<-row+i
                test_names<-row.names(dist1[order(dist1[,start]),])[1:nn.number]
                train_set.mw<-train_set[test_names,]
                train_env.mw<-train_env[test_names,]
                test_set.mw<-as.matrix(test_set[i,])
                
                
                if (rplot==TRUE && mw.type!="cca")
                {
                    points(test_scores[i,1],test_scores[i,2],col="red",pch=19,cex=2)     
                    points(train_scores[,1][row.names(train_scores)%in%test_names],train_scores[,2][row.names(train_scores)%in%test_names],col="blue",pch=19)
                
                par(fig=c(0,1,0.6,1),new=TRUE)
                par(mar=c(0,5,0,0))
                plot(train_scores[,1],train_scores[,2],type="n",axes=FALSE,xlab="",ylab="")
                legend("topleft",c("training set","MW-training set", "core samples", "finish core samples"),pt.cex=1,cex=1,pch=c(21,19,19,19),col=c(1,"blue","green","orange"),box.lty=0,ncol=2)
                }
                
                
                
                if (method=="wapls")
                    {
                        mwstat<-wapls(train_set.mw,train_env.mw,test_set,d.plot="FALSE",diagno="FALSE",comp=ncomp,val=val,run=run,out="FALSE",drop.non.sig=drop.non.sig,min.occ=min.occ,scale=scale)
                        wapls_comp<-best.wapls(mwstat$performance,nro=ncomp,rmsep.incl=rmsep.incl)
                        wapls_comp.mw[k,]<-wapls_comp[ncomp,]
                        mw_env[k]<-mwstat$reconstruction_core.samples[i,wapls_comp[ncomp,9]]
                        if (val%in%c("boot","loo"))
                        { 
                            mw_env.val[k]<-mwstat$"mean(reconstruction_core.samples).val"[i,wapls_comp[ncomp,9]]    
                            mw_env.sd[k]<-mwstat$"sd(reconstruction_core.samples).val"[i,wapls_comp[ncomp,9]]      
                        }
                   
                    }
                
                if (method=="wa")
                    {
                        mwstat<-wa(train_set.mw,train_env.mw,test_set,d.plot="FALSE",diagno="FALSE",val=val,run=run,out="FALSE",drop.non.sig=drop.non.sig,min.occ=min.occ,scale=scale)
                        wa_per<-c(mwstat$performance,1,1)
                        wapls_comp.mw[k,]<-wa_per
                        mw_env[k]<-mwstat$reconstruction_core.samples[i]
                        if (val%in%c("boot","loo"))
                            { 
                                mw_env.val[k]<-mwstat$"mean(reconstruction_core.samples).val"[i]   
                                mw_env.sd[k]<-mwstat$"sd(reconstruction_core.samples).val"[i]     
                            }
                    }
          
                if (method=="pom")
                    {
                        mwstat<-pom(train_set.mw,train_env.mw,test_set,d.plot="FALSE",val="loo",out="FALSE",scale=TRUE)
                        wa_per<-c(mwstat$performance,1,1)
                        wapls_comp.mw[k,]<-wa_per
                        mw_env[k]<-mwstat[[5]][i]
                        if (val%in%c("boot","loo"))
                            { 
                                mw_env.val[k]<-mwstat[[6]][i]   
                                mw_env.sd[k]<-mwstat[[7]][i]     
                            }
                    }
          
          
          
          
          
          
           }
        mw_env
        comp.r<-wapls_comp.mw[,9]
        wapls_comp.mw[,9]<-c(1:number.mw_ds)
        end_set<-best.wapls(wapls_comp.mw,nro=(number.mw_ds),rmsep.incl=rmsep.incl)
        end_set
        mw_env.end[i]<-mw_env[end_set[number.mw_ds,9]]
        if (val%in%c("boot","loo"))
            {
                mw_env.end.val[i]<-mw_env.val[end_set[number.mw_ds,9]]
                mw_env.end.sd[i]<-mw_env.sd[end_set[number.mw_ds,9]]
            }
        end_set_v<-as.vector(end_set[number.mw_ds,])
        mw_envstat.end[i,]<-c(round(end_set_v[1:8],4),round(comp.r[end_set[number.mw_ds,9]],0),round(mwsize[end_set[number.mw_ds,9]],0))
    
    
    
    
    }


res<-list(mw_envstat.end,mw_env.end,mw_env.end.val,mw_env.end.sd)
names(res)[[1]]<-"sample.performance"
names(res)[[2]]<-"reconstruction"
names(res)[[3]]<-"mean(reconstruction).val"
names(res)[[4]]<-"sd(reconstruction).val"
inf_test<-mw_env.end
x<-c(1:length(inf_test))
res
}