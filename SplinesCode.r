require(gss)
require(pracma)


##permTestNew function computes the 1000 different permutations for class

permTestNew<-function(dat,t=0,B=1000){
  counts = dat$resp;
  time   = dat$posi;
  statusP= dat$cl;
  mat = array(dat$cl,dim=c(length(data$cl),B))
  tmp = strsplit(rownames(dat),split="")
  #tmp2 = sapply(tmp,function(i){
    #x = i[[1]][1]
      #for(j in 2:4){ x=paste(x,i[j],sep="")}
    #return(x)})
  tmp2=dat$inde
  
  k = match(unique(tmp2),tmp2)
  cl  = statusP[k+3];
  
  for(l in 1:B){
    clp = sample(cl,length(cl));

    for(j in 1:(length(k)-1)){
      statusP[k[j]:k[j+1]]=clp[j]
    }
    statusP[k[j+1]:length(statusP)]=clp[j+1]
    
    mat[,l] = statusP;
  }
  dat$statusp = mat;
  
  return(dat)
}


####SSRegionFinder uses SSANOVA to find regions of difference. 
###Function takes in data with response, class, Individual ID, position, and perm # for to perform


SSRegionFinder=function(data, Y, cl, indID, position, Perm){

require(gss)
require(pracma)

#resp=data[,Y]
#cla=factor(data[,cl])
#inde=data[,indID]
#posi=data[,position]
dat=data.frame(resp=data[,Y], cla=factor(data[,cl]), inde=data[,indID], posi=data[,position])
mod=ssanova(resp~posi*cla,partial=~inde,data=dat)
x_seq=seq(min(dat$posi),max(dat$posi),by=1)
predP=predict(mod,data.frame(posi=x_seq, cla=factor(1)), se=T,  inc=c("cla","posi:cla"))

Region_pos=which((2*predP$fit+(1.645*2*predP$se))<0)

if(length(Region_pos)>0){

if(sum(diff(Region_pos)>1)<1){
	k_pos=list()
	k_pos[[1]]=Region_pos
}else{
	k_pos=list()
	useInd=which(diff(Region_pos)>1)
	i=1
	while (i<=length(useInd))
		{
			k_pos[[i]]=Region_pos[1:useInd[i]]
			Region_pos=Region_pos[-(1:useInd[i])]	
			i=i+1	
		}
	i=length(k_pos)+1
	k_pos[[i]]=Region_pos
	}
}else{
	k_pos=c()
	}

Region_neg=testpl=which((2*predP$fit-(1.645*2*predP$se))>0)
if(length(Region_neg)>0){

if(sum(diff(Region_neg)<1)){
	k_neg=list()
	k_neg[[1]]=Region_neg
}else{
	k_neg=list()
	useInd=which(diff(Region_neg)>1)
	i=1
	while (i<=length(useInd))
		{
			k_neg[[i]]=Region_neg[1:useInd[i]]
			Region_neg=Region_neg[-(1:useInd[i])]	
			i=i+1	
		}
	i=length(k_neg)+1
	k_neg[[i]]=Region_neg
	}
}else{
	k_neg=c()
	}

k=c(k_pos,k_neg)

#return(k)
#}

#if (k>0){

Result=matrix(0, nrow=length(k), ncol=4)
colnames(Result)=c("Start", "End", "Area", "P-Value")

for (i in 1:nrow(Result)){
Result[i,1]=min(k[[i]])
Result[i,2]=max(k[[i]])
Result[i,3]=trapz(x=x_seq[min(k[[i]]):max(k[[i]])], y=abs(2*predP$fit[min(k[[i]]):max(k[[i]])]))
}

new_dat=permTestNew(dat, t=0, B=Perm)
areaPSp=matrix(0,Perm,nrow(Result))
dat2=dat

#fine up to here
#return(Result)
#}


for (i in 1:Perm){
        
        dat2$cla=factor(as.numeric(new_dat$statusp[,i])-1)
        mod=ssanova(resp~posi*cla, partial=~inde, data=dat2)
        predPl=predict(mod, data.frame(posi=x_seq, cla=factor(1)), se=T, inc=c("posi","posi:cla"))
        	for (j in 1:nrow(Result)){
							areaPSp[i, j]=trapz(x=x_seq[Result[j,1]:Result[j,2]],  y=abs(2*predPl$fit[Result[j,1]:Result[j,2]]))
						}				        
        #show(i)
}

for (i in 1:nrow(Result)){
	Result[i,4]=1-mean(Result[i,3]<areaPSp[,i])
}

#Here its fine
if (length(k)==0){
return("No Region Found")}else{
	return(Result)}

}


#Example: SSRegionFinder(dat, 1, 2, 3, 4, 20)
