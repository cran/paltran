`palplot` <-
function(...,age=NA,error.plot=FALSE,error=0,ptype="b",error.ptype="b",p.col="black",p.lty=1,p.xlab=NULL,cex.main=1,cex.axis=1,dis=0.15,p.main=NULL,y.lab="sample",trend=0,span=0.1)
{
    par(mfrow=c(1,1))
    split.screen(c(1,1))
    data<-list(...)
    data_ges<-NA
    data_l<-NA
    for (i in (seq_len(length(data))))
        {   if (is.data.frame(data[[i]])==TRUE)
                data_l<-as.list(data[[i]])
            if (is.matrix(data[[i]])==TRUE)
                data_l<-as.list(as.data.frame(data[[i]]))
            if (is.vector(data[[i]])==TRUE)
                data_l<-as.data.frame(data[[i]])
            if (i==1)
                data_ges<-data_l
            if (i>1)
                data_ges<-c(data_ges,data_l)
        }

    if (is.na(age[1]))
        y<-seq_len(length(data_ges[[1]]))
        else
            y<-age

    len<-length(data_ges)
    dist<-(1-0.08)/(len)
    par(mfrow=c(1,1))
    par(fig=c(0,0.01,0,1))
    if (length(ptype)<len)
        sym<-rep(ptype,len)
        else
            sym<-ptype
    if (length(p.col)<len)
        p.col1<-rep(p.col,len)
        else
            p.col1<-p.col
    if (length(p.xlab)<len)
        p.xlab1<-rep(p.xlab,len)    
        else
            p.xlab1<-p.xlab
if (length(trend)<len)
        trend<-rep(trend,len)
        else
            trend<-trend

error<-as.matrix(error)

    nr<-NA
    for (k in seq_len(len))
    {
        par(mgp=c(2,0.5,0))
        start<-(k-1)*dist+0.07
        end<-(k*dist)+0.07
        par(fig=c(start,end,0,1),new=TRUE)
        if (k ==1)
            par(mar=c(4,.2,3,dis))
        if (k>1)
            par(mar=c(4,dis,3,dis))
        x<-as.vector(data_ges[[k]])
        spec<-x
        leange<-length(x)
        
        
        
        
        if (sym[k]=="h")
            plot(x,-y,las=1,ylim=c(min(-y),max(-y)),yaxt="n",type="n",xlim=c(0,max(x)+0.5*max(x)),xlab=p.xlab1[k],main=p.main[k],cex.main=cex.main,cex.axis=cex.axis)
            
        if (sym[k]=="b")
         {
                plot(x,-y,las=1,ylim=c(min(-y),max(-y)),yaxt="n",type="b",col=p.col1[k],lty=p.lty,xlim=c(0,max(x)+0.5*max(x)),xlab=p.xlab1[k],main=p.main[k],cex.main=cex.main,cex.axis=cex.axis)
                if (error.plot=="TRUE")
                {
                    if (error.ptype=="b")
                        if (sum(abs(error[,k]))!=0)
                        {
                          for (i in seq_len(length(x)))
                        {
                        arrows(x[i],-y[i],x[i]+error[i,k],-y[i],angle=90,length = 0.05)
                        arrows(x[i],-y[i],x[i]-error[i,k],-y[i],angle=90,length = 0.05)
                        }
                    }
                    if (error.ptype=="p")
                        if (sum(abs(error[,k]))!=0)
                        {
                        xp1<-x
                        xp2<-x-error[,k]
                        xp2.1<-x+error[,k]
                        yp1<--y
                        yp2<--y
                        polygon(c(xp1[1],xp2,xp1[order(yp1)]),c(yp1[1],yp2,yp1[order(yp1)]),col="lightgray",border="lightgray")
                        polygon(c(xp1[1],xp2.1,xp1[order(yp1)]),c(yp1[1],yp2,yp1[order(yp1)]),col="lightgray",border="lightgray")
                        points(x,-y,col=p.col1[k],lty=p.lty,type="b",pch=19,col=p.col)
                        }
                }
        
            }
        
        if (sym[k]=="l")
            {
                plot(x,-y,las=1,ylim=c(min(-y),max(-y)),yaxt="n",type="l",col=p.col1[k],lty=p.lty,xlim=c(0,max(x)+0.5*max(x)),xlab=p.xlab1[k],main=p.main[k],cex.main=cex.main,cex.axis=cex.axis)
                if (error.plot=="TRUE")
                {
                    if (error.ptype=="b")
                        if (sum(abs(error[,k]))!=0)
                        {
                            for (i in seq_len(length(x)))
                        {
                        arrows(x[i],-y[i],x[i]+error[i,k],-y[i],angle=90,length = 0.05)
                        arrows(x[i],-y[i],x[i]-error[i,k],-y[i],angle=90,length = 0.05)
                        }
                        }
                    if (error.ptype=="p")
                        
                        {
                        xp1<-x
                        xp2<-x-error[,k]
                        xp2.1<-x+error[,k]
                        yp1<--y
                        yp2<--y
                        polygon(c(xp1[1],xp2,xp1[order(yp1)]),c(yp1[1],yp2,yp1[order(yp1)]),col="lightgray",border="lightgray")
                        polygon(c(xp1[1],xp2.1,xp1[order(yp1)]),c(yp1[1],yp2,yp1[order(yp1)]),col="lightgray",border="lightgray")
                        points(x,-y,col=p.col1[k],lty=p.lty,type="l",pch=19)
                        }
                }
        
            }
        
        if (sym[k]=="p")
            {
            plot(x,-y,las=1,ylim=c(min(-y),max(-y)),yaxt="n",type="p",col=p.col1[k],lty=p.lty,xlim=c(0,max(x)+0.5*max(x)),xlab=p.xlab1[k],main=p.main[k],cex.main=cex.main,cex.axis=cex.axis)
            if (error.plot=="TRUE")
                {
                    if (error.ptype=="b")
                        if (sum(abs(error[,k]))!=0)
                        {
                        for (i in seq_len(length(x)))
                        {
                        arrows(x[i],-y[i],x[i]+error[i,k],-y[i],angle=90,length = 0.05)
                        arrows(x[i],-y[i],x[i]-error[i,k],-y[i],angle=90,length = 0.05)
                        }
                    }
                    if (error.ptype=="p")
                        
                        {
                        xp1<-x
                        xp2<-x-error[,k]
                        xp2.1<-x+error[,k]
                        yp1<--y
                        yp2<--y
                        polygon(c(xp1[1],xp2,xp1[order(yp1)]),c(yp1[1],yp2,yp1[order(yp1)]),col="lightgray",border="lightgray")
                        polygon(c(xp1[1],xp2.1,xp1[order(yp1)]),c(yp1[1],yp2,yp1[order(yp1)]),col="lightgray",border="lightgray")
                        points(x,-y,col=p.col1[k],lty=p.lty,type="l",pch=19)
                        }
                }
            }
        if (trend[k]==1)
            {
                fit.loe<-loess(x~y,span=span)
                ifelse(fitted(fit.loe)<0,0,fitted(fit.loe))
                points(fitted(fit.loe),-y,type="l",col="red",lwd=2)
            }
        
        if (k==1)
        {
            nr<-NA
            anz<-round(length(x)/10,0)
            for (i in 1:(anz+1))
                if (i==1)
                    nr[i]<-1
                else
                    nr[i]<-10*(i-1)
            axis(2,at=-y[nr],labels=round(y[nr],0),las=1,cex.axis=cex.axis)
            mtext(side=2,text=y.lab,line=3,font=2)
        }
        x2<-spec
        y1<-(-y)#(seq(1:leange))
        x1<-rep(0,leange)
        y2<-(-y)#(seq(1:leange))
        if (sym[k]=="h")
        {
            for (i in 1:leange)
                lines(cbind(x1,x2)[i,],cbind(y1,y2)[i,],lwd=3,col=p.col1[k])
        }


    }
    close.screen (all=TRUE)
}

