
library(gss)
library(pracma)

load(SampleData.rDa)

#Make Data, Run original data and produce fit and standard errors

SSfit <- function(data, response, group, series, id)
	{
		dat=data.frame(response=data[,response], group=factor(data[,group]), series=data[,series], id=data[,id])
		mod=ssanova(response ~ series * group, data=dat)
		full_x = seq(min(dat$series), max(dat$series), by=1)
		fit = predict(mod, data.frame(series=full_x, group=factor(1)), inc=c("group", "series:group"), se=TRUE)
		return(list(data=dat, fit=fit$fit, se=fit$se, x_seq=full_x))
	}


#Permutation matrix (Input: Ready made data)
SSmakePerm=function(Data_Ready)
	{
		dat2=data.frame(group=factor(Data_Ready$group), id=Data_Ready$id)
		count=table(dat2$id)
		samp=unique(dat2)[,1]
		permMat=list()
		
		for (j in 1:1000){
			new_samp=sample(samp, replace=FALSE)
			permMat[[j]]=unlist(sapply(1:length(new_samp), function(i){rep(new_samp[i], count[i])}))
		}
		return(permMat)	
	}

#Analysis of all permutations
SSPermAnalysis=function(Data_Ready, PermList, ResultMatrix, x_seq){
	result_perm=matrix(0, 1000, nrow(ResultMatrix))
	x_perm=x_seq
	dat_perm=Data_Ready

		for (j in 1:1000){
			dat_perm$group=PermList[[j]]
			mod_perm=ssanova(response ~ series * group, data=dat_perm)
			fit_perm = cbind(x_perm, abs(2*predict(mod_perm, data.frame(series=x_perm, group=factor(1)), inc=c("group", "series:group"), se=TRUE)$fit))
				for (i in 1:nrow(ResultMatrix)){
					area_perm=fit_perm[which(fit_perm[,1]==ResultMatrix[i,1]) : which(fit_perm[,1]==ResultMatrix[i, 2]), ]
					result_perm[j, i]=trapz(x=area_perm[,1], y=area_perm[,2])
				}
			if (j%%100==0){
				show(j)
			}	
		}
	return(result_perm)
}


#Check if any region where curve is above 0
#input (fit, se, sequence, above or below line 0)
SSRegions_candidate=function(fit, stand.error, x_seq, above_zero=TRUE){
	if (above_zero){
	diff_position = which((2*fit-(1.96*2*stand.error))>0)
	}else{
		diff_position = which((2*fit+(1.96*2*stand.error))<0)
	}
	if (length(diff_position)>0){
		reg_index=which(diff(diff_position)!=1)
		reg_position=matrix(0, (length(reg_index)+1), 4)
			if (length(reg_index)==0){
				reg_position[1,1]=x_seq[diff_position[1]]
				reg_position[1,2]=x_seq[tail(diff_position, n=1)]
			}else{
				i=1
				while(length(reg_position)!=0 & length(reg_index)!=0){
					reg_position[i, 1]=x_seq[diff_position[1]]
					reg_position[i,2]=x_seq[diff_position[reg_index[1]]]
					diff_position=diff_position[-c(1:reg_index[1])]
					reg_index=reg_index[-1]
					i=i+1
				}
			reg_position[i, 1]=x_seq[diff_position[1]]
			reg_position[i, 2]=x_seq[tail(diff_position, n=1)]
			}
		}else{
			reg_position=NULL	
	}	
	return(reg_position)	
}


SSRegion_Finder=function(data, response, group, series, id ){
	set.seed(123)
	dat_original=SSfit(data=data, response=response, group=group, series=series, id=id)
	#perm_mat=SSmakePerm(dat_original$data)
	index_pos=SSRegions_candidate(fit=dat_original$fit, stand.error=dat_original$se, x_seq=dat_original$x_seq, above_zero=TRUE)
	index_neg=SSRegions_candidate(fit=dat_original$fit, stand.error=dat_original$se, x_seq=dat_original$x_seq, above_zero=FALSE)
	index_all=rbind(index_pos, index_neg)
	
	if (nrow(index_all)>0){
		colnames(index_all)=c("Region Start", "Region End", "Area", "P-Value")
		predict_area=cbind(dat_original$x_seq , abs(2*dat_original$fit))
		perm_list=SSmakePerm(dat_original$data)
		PermResult=SSPermAnalysis(Data_Ready= dat_original$data, PermList=perm_list, ResultMatrix=index_all, x_seq=dat_original$x_seq )
		for (i in 1:nrow(index_all)){
			area_orig=predict_area[which(predict_area[,1]==index_all[i,1]):which(predict_area[,1]==index_all[i, 2]), ]
			actual_area=trapz(x=area_orig[,1], y=area_orig[,2])
			index_all[i, 3]=actual_area
			index_all[i, 4]=1-(length(which(actual_area>PermResult[,i]))/1000)
		}
	
	}else{
		return("No Regions Found")}
	
	return(list(Result=index_all, data.frame=dat_original$data, fit=2*dat_original$fit, se=2*dat_original$se))
}


# Test function RUN THIS#####

res=SSRegion_Finder(data=data, response=3, group=2, series=1, id=4)
