`best.wapls` <-
function(error,nro=4,rmsep.incl=TRUE)
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
                error[,9]<-c(1:nro)
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

