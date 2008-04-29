`mwtraf` <-
function(...,method="wapls",rplot="TRUE",mwsize=c(40,60,80,100,120,140,160),mw_type="dca",ncomp=4,rmsep.incl=TRUE,env.trans=FALSE,spec.trans=FALSE)

{
data<-list(...)
train_set<-data[[1]]
train_env<-data[[2]]


test_set<-data[[3]]


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

best.wapls<-function(error,nro=4)
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
wapls_comp.mw<-matrix(nrow=(number.mw_ds),ncol=10)
nam1<-row.names(test_set)
mw_envstat.end<-matrix(NA,nrow=row1,ncol=10,dimnames=list(nam1,c("RMSE","Rs","Ave_Bias","Max_Bias","Rs.X","Ave_Bias.X","Max_Bias.X","RMSEP","comp","w.size")))
mw_env.end<-rep(0,row1)
new<-merge(train_set,test_set,all=TRUE,sort=FALSE)                      # creating test and training_set with the same rows
new<-as.data.frame(new)
row.names(new)<-c(row.names(train_set),row.names(test_set))
new<-new[,(1:col1)]
train_set<-new[(1:row),]
train_set<-apply(train_set,2,function(x){ifelse(is.na(x),0,x)})
test_set<-new[((row+1):dim(new)[1]),]
test_set<-apply(test_set,2,function(x){ifelse(is.na(x),0,x)})
wapls_comp<-matrix(0,ncol=10,nrow=ncomp)


######################################################################################

library(vegan)
if (mw_type=="dca")
{
train.dca<-decorana(downweight(train_set)) 
train_scores<-scores(train.dca)
test.dca<-predict(train.dca,newdata=test_set,type="sites")  
test_scores<-scores(test.dca)
test_train.co<-rbind(train_scores,test_scores)
dist<-vegdist(test_train.co,method="euclidean" )
dist<-as.matrix(dist)
dist<-as.data.frame(dist)

}
row<-dim(train_set)[1]
dist1<-dist[1:row,]


for(i in 1:row1)
{
for (k in 1:(number.mw_ds))
{
start<-row+i
nn.number<-mwsize[k]
test_names<-row.names(dist1[order(dist1[,start]),])[1:nn.number]

train_set.mw<-train_set[test_names,]
train_env.mw<-train_env[test_names,]
test_set.mw<-as.matrix(test_set[i,])

mwstat<-wapls(train_set.mw,train_env.mw,test_set,d.plot=FALSE,diagno=FALSE,n_comp=ncomp)
wapls_comp<-best.wapls(mwstat$performance,nro=ncomp)
wapls_comp.mw[k,]<-wapls_comp[ncomp,]

mw_env[k]<-mwstat$reconstruction[i,wapls_comp[ncomp,9]]
}
mw_env
wapls_comp.mw
end_set<-best.wapls(wapls_comp.mw,nro=(number.mw_ds))
end_set
mw_env.end[i]<-mw_env[end_set[number.mw_ds,9]]
end_set_v<-as.vector(end_set[number.mw_ds,])
mw_envstat.end[i,]<-c(round(end_set_v[1:9],4),round(mwsize[end_set[number.mw_ds,9]],0))
}
res<-list(mw_envstat.end,mw_env.end)
names(res)[[1]]<-"sample.performance"
names(res)[[2]]<-"reconstruction"
        inf_test<-mw_env.end
        x<-c(1:length(inf_test))
        plot(-x~inf_test,axes=FALSE,type="l",col="darkred",ylab="sample ",main="MW - reconstruction",xlab="environmental parameter",xlim=c(0,max(inf_test)+1))#
        box()
        axis(1)
        points(inf_test,-x,type="p",pch=19,cex=0.5)
        segments(inf_test,-x,(inf_test-mw_envstat.end[,8]),-x,col="red")
        segments(inf_test,-x,(inf_test+mw_envstat.end[,8]),-x,col="red")
res
}
