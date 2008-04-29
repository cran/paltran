rand.test <-
function(error1,error2,ran.nb=1000,seed=1)

{
result<-matrix(ncol=3,nrow=1,dimnames = list(c(""),c("MSEP1","MSEP2","p")))
e1<-error1
e2<-error2
a<-length(e1)
b<-length(e2)
if (a!=b)
   {
   print(paste("unequal dimension of the errors! using only the first",a,"to calculate p"))
   e2<-e2[1:a]
   }
d<-e1^2-e2^2
MSEP1<-sum(e1^2)/a
MSEP2<-sum(e2^2)/a

T.ges<-mean(d)
set.seed(seed)
sum(d)/a

Ti<-NA
s<-0
for (j in 1:ran.nb)
    {
    sig<-NA
    sig<-sample(c(1,-1),a,replace=TRUE)
    d1<-ifelse(sig>0,abs(d),-abs(d))
    Ti[j]<-mean(d1)
    if (abs(Ti[j])>=abs(T.ges))
       s<-s+1
    }
p<-(s)/(ran.nb+1)
result[1,]<-round(c(MSEP1,MSEP2,p),3)

result
}

