require(gss)
require(pracma)


###data Frame with 3 columns, 1 for location in the genomoe, 1 methylation level and 1 binary for tumor or not (1, 0)


# Splines#########
mod=ssanova(methy~loca*case,data=dat3)

x=seq(min(dat3$loca),max(dat3$loca),by=1)

pred1=predict(mod,data.frame(loca=x,case=factor(1)),se=T)$fit
pred2=predict(mod,data.frame(loca=x,case=factor(2)),se=T)$fit


predP=predict(mod,data.frame(loca=x,case=factor(2)),se=T, inc=c("case","loca:case"))

#Find where difference function is less than 0
testpl=which((2*predP$fit+(1.645*2*predP$se))<0)

k=which(diff(testpl)!=1)
k1=testpl[1:k[1]]
k2=testpl[(k[1]+1):(k[2])]

#find where difference function is greater than 0
testmin=which((2*predP$fit-(1.645*2*predP$se))>0)
diff(testmin)

#Calculate real areas
areaA=trapz(x=x[k1],y=abs(2*predP$fit[k1]))
areaB=trapz(x=x[testmin],y=abs(2*predP$fit[testmin]))
areaC=trapz(x=x[k2],y=abs(2*predP$fit[k2]))


#1000 permutations and calculate areas
head(dat3)
dat4=dat3
areaPSp=matrix(0,1000,3)

for (i in 1:1000){
	
	dat4$case=rep(sample(tumor),each=22)
	mod=ssanova(methy~loca*case,data=dat3)
	predPl=predict(mod,data.frame(loca=x,case=factor(2)),se=T,inc=c("case","loca:case"))
	areaPSp[i,1]=trapz(x=x[k1],y=abs(2*predPl$fit[k1]))
	areaPSp[i,2]=trapz(x=x[testmin],y=abs(2*predPl$fit[testmin]))
	areaPSp[i,3]=trapz(x=x[k2],y=abs(2*predPl$fit[k2]))
	show(i)
	
}

#Compares Calc areas with permutations .95%
