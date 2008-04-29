`palplot` <-
function(...,col="red",trend=TRUE,trendcol="gray90",lty=1,lwd=1,k=FALSE,polyg=TRUE,colpoly="gray90",trans=FALSE,meth="loess",span=0.2,conf.t=TRUE,xlab="inferred env. factor", ylab="sample",main="reconstruction",cex=1.5,pch=20,type="p")

{
data<-list(...)
env<-data[[1]]

if (length(data)==2)
        age<-data[[2]]  
          
if (length(data)<=1)
        age<-c(1:length(env))
        

ord<-order(age)
age<-age[ord]
env<-env[ord]


if (trans=="log10")      
    env<-10^env
if (trans=="log")
    env<-exp^env
if (trans=="sqrt")
    env<-env^2


plot(env,-age,pch=pch,ylab=ylab,xlab=xlab,main=main,ylim=c(min(-age),max(-age)))
        
conf.up<-NA
conf.un<-NA
sum<-NA
m.trend<-NA

if (trend==TRUE)
    {
        if(meth=="gam")
            {
                library(mgcv)
                if (k==FALSE)
                    k<-round(length(env)/10,0)
                fit.gam<-gam(env~s(age,k=k))
                
                if (trend==TRUE)
                    {
                        fit.pre<-predict(fit.gam,se=TRUE)
                        points(fit.pre$fit+1.96*fit.pre$se.fit,-age,type="l",col=trendcol)
                        points(fit.pre$fit-1.96*fit.pre$se.fit,-age,type="l",col=trendcol)
                        conf.up<-fit.pre$fit+1.96*fit.pre$se.fit
                        conf.un<-fit.pre$fit-1.96*fit.pre$se.fit
                        if (polyg==TRUE)
                            {
                                n<-length(env)
                                x<-c(-age,-age[n:1])
                                y<-c(conf.un,conf.up[n:1])
                                polygon(y,x,col=colpoly,border=NA)
                            }
                    }  
                points(fitted(fit.gam),-age,type="l",col=col,lty=lty,lwd=lwd)
                m.trend<-fitted(fit.gam)
                sum<-summary(fit.gam)
            }
        
        if (meth=="loess")
            {
                fit.loe<-loess(env~age,span=span)
                points(fitted(fit.loe),-age,type="l",col=col)
                if (trend==TRUE)
                    {
                        fit.pre<-predict(fit.loe,se=TRUE)
                        points(fit.pre$fit+1.96*fit.pre$se.fit,-age,type="l",col=trendcol)
                        points(fit.pre$fit-1.96*fit.pre$se.fit,-age,type="l",col=trendcol)
                        conf.up<-fit.pre$fit+1.96*fit.pre$se.fit
                        conf.un<-fit.pre$fit-1.96*fit.pre$se.fit
                        if (polyg==TRUE)
                            {
                                n<-length(env)
                                x<-c(-age,-age[n:1])
                                y<-c(conf.un,conf.up[n:1])
                                polygon(y,x,col=colpoly,border=NA)
                            }
                    }       
                points(fitted(fit.loe),-age,type="l",col=col,lty=lty,lwd=lwd)
                m.trend<-fitted(fit.loe)
                sum<-summary(fit.loe)
            }
        if (meth=="scatter")
        {}
            
        
    }


points(env,-age,pch=pch,ylab=ylab,cex=cex)
box()
results<-list(env,m.trend,conf.up,conf.un,sum)
names(results)[[1]]<-"env.inferred"
names(results)[[2]]<-"env.trend"
names(results)[[3]]<-"trend.conf.up"
names(results)[[4]]<-"trend.conf.un"
names(results)[[5]]<-"summary.trend"
results
}

