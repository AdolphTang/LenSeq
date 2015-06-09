################################################################################################################
binomTest <- function(y1_bino, y2_bino, n1=sum(y1_bino), n2=sum(y2_bino), p=n1/(n1+n2))
{
	if(length(y1_bino) != length(y2_bino)) stop("y1_bino and y2_bino must have same length")
	if(any(is.na(y1_bino)) || any(is.na(y2_bino))) stop("missing values not allowed")
	y1_bino <- round(y1_bino)
	y2_bino <- round(y2_bino)
	if(any(y1_bino<0) || any(y2_bino<0)) stop("y1_bino and y2_bino must be non-negative")
	if(p<=0 || p>=1) stop("p must be between 0 and 1")
	size <- y1_bino+y2_bino
	p.value <- rep.int(1,length(y1_bino))
	if(p==0.5) {
		i <- (size>0)
		if(any(i)) {
			y1_bino <- pmin(y1_bino[i],y2_bino[i])
			size <- size[i]
			p.value[i] <- pmin(2*pbinom(y1_bino,size=size,prob=0.5),1)
		}
		return(p.value)
	}
	if(any(big <- size>10000)) {
		ibig <- which(big)
		for (i in ibig) p.value[i] <- chisq.test(matrix(c(y1_bino[i],y2_bino[i],n1-y1_bino[i],n2-y2_bino[i]),2,2))$p.value
	}
	size0 <- size[size>0 & !big]
	if(length(size0)) for (isize in unique(size0)) {
		i <- (size==isize)
		d <- dbinom(0:isize,prob=p,size=isize)
		o <- order(d)
		cumsump <- cumsum(d[o])[order(o)]
		p.value[i] <- cumsump[y1_bino[i]+1]
	}
	p.value
}

exactTest <- function(y1,y2,dispersionA=0,dispersionB=0,big.count=6000)
{
#	Convert matrices to vectors
	ntags <- NROW(y1)
	n1 <- NCOL(y1)
	n2 <- NCOL(y2)
	if(n1>1) s1 <- round(rowSums(y1)) else s1 <- round(y1)
	if(n2>1) s2 <- round(rowSums(y2)) else s2 <- round(y2)
	s <- s1+s2
	mu <- s/(n1+n2)
	mu1 <- n1*mu
	mu2 <- n2*mu

	pvals <- rep(1,ntags)
	names(pvals) <- names(y1)

#	Poisson case
	poisA <- dispersionA<=0
	poisB <- dispersionB<=0
	pois = poisA | poisB
#	BINOMTEST DOESN'T USE EQUAL TAILED REJECTION REGION
	if(any(pois)) pvals[pois] <- binomTest(s1[pois],s2[pois])


	p.bot <- size1 <- size2 <- rep(0,ntags)
	left <- s1<mu1 & !pois# & !big        
	if(any(left)) {
		p.bot[left] <- dnbinom(s[left],size=(n1+n2)/dispersionB[left],mu=s[left])   # dispersionA or B ?????????????????
		size1[left] <- n1/dispersionA[left]
		size2[left] <- n2/dispersionB[left]
		for (g in which(left)) {
			x <- 0:s1[g]
			p.top <- dnbinom(x,size=size1[g],mu=mu1[g]) * dnbinom(s[g]-x,size=size2[g],mu=mu2[g])
			pvals[g] <- 2*sum(p.top)
		}
		pvals[left] <- pvals[left]/p.bot[left]
	}
	right <- s1>mu1 & !pois# & !big
	if(any(right)) {
		p.bot[right] <- dnbinom(s[right],size=(n1+n2)/dispersionA[right],mu=s[right])
		size1[right] <- n1/dispersionA[right]
		size2[right] <- n2/dispersionB[right]
		for (g in which(right)) {
			x <- s1[g]:s[g]
			p.top <- dnbinom(x,size=size1[g],mu=mu1[g]) * dnbinom(s[g]-x,size=size2[g],mu=mu2[g])
			pvals[g] <- 2*sum(p.top)
		}
		pvals[right] <- pvals[right]/p.bot[right]
	}
	pmin(pvals,1)
}

estSizeFactors <- function( counts, locfunc = median )
{
   loggeomeans <- rowMeans( log(counts) )
   apply( counts, 2, function(cnts)
      exp( locfunc( ( log(cnts) - loggeomeans )[ is.finite(loggeomeans) ] ) ) )
}

normalized <- function(counts , sizeFac)
{
    t(t(counts)/sizeFac)
}

############################################################################################
# Inputs
############################################################################################
data = read.table("read_count.txt",header=TRUE)
geneLength = read.table("geneLen_from_gtf.txt",header=FALSE)

############################################################################################
# Normalization
############################################################################################
#gene length match
match_index = match(geneName,geneLength[,1])
match_index = na.omit(match_index)
geneLength = geneLength[match_index,]
length_mean = mean(geneLength[,2])

#this 0.5parameter can be adjust
geneLength[,2] = 0.5*geneLength[,2]+0.5*length_mean
length_mean = mean(geneLength[,2])
################################################
match_index = match(geneLength[,1],rownames(data))
match_index = na.omit(match_index)

geneLength_data = rep(NA,nrow(data))
geneLength_data[match_index] = geneLength[,2]
geneLength_data[which(is.na(geneLength_data))] = length_mean

library(matrixStats)

#size factor
sf=estSizeFactors(data)

sizeB = sf[1:(ncol(data)/2)]
sizeU = sf[(ncol(data)/2+1):ncol(data)]

data_norm = normalized(data,sf)

mean_B = rowMeans(data_norm[,1:(ncol(data)/2)])
var_B = rowVars(data_norm[,1:(ncol(data)/2)])
mean_U = rowMeans(data_norm[,(ncol(data)/2+1):ncol(data)])
var_U = rowVars(data_norm[,(ncol(data)/2+1):ncol(data)])
# calculate raw dispersion
ximB <- mean( 1/sizeB )
disp_B <- ( var_B - ximB * mean_B ) / mean_B^2
ximU <- mean( 1/sizeU )
disp_U <- ( var_U - ximU * mean_U ) / mean_U^2

mean_B_tmp = mean_B
mean_U_tmp = mean_U
disp_B_tmp = disp_B
disp_U_tmp = disp_U
mean_B = mean_B_tmp[mean_B_tmp>0]
disp_B = disp_B[mean_B_tmp>0]
disp_U = disp_U[mean_U_tmp>0]
mean_U = mean_U_tmp[mean_U_tmp>0]

library(edgeR)
data_div_norm = (data_norm*length_mean)/geneLength_data

input.mean = as.matrix(data_norm[mean_B_tmp>0,1:(ncol(data)/2)])
output.mean = as.matrix(data_div_norm[mean_B_tmp>0,1:(ncol(data)/2)])

disp_B[which(disp_B<0)] = 0.0000008
test = q2qnbinom(input.mean,input.mean=input.mean,output.mean=output.mean,dispersion=disp_B)
data_div_norm[mean_B_tmp>0,1:(ncol(data)/2)] = test
input.mean = data_norm[mean_U_tmp>0,(ncol(data)/2+1):ncol(data)]
output.mean = data_div_norm[mean_U_tmp>0,(ncol(data)/2+1):ncol(data)]
disp_U[which(disp_U<0)] = 0.0000008
test = q2qnbinom(input.mean,input.mean=input.mean,output.mean=output.mean,dispersion=disp_U)
data_div_norm[mean_U_tmp>0,(ncol(data)/2+1):ncol(data)] = test
rm(test)
mean_B = rowMeans(data_div_norm[,1:(ncol(data)/2)])
mean_U = rowMeans(data_div_norm[,(ncol(data)/2+1):ncol(data)])
mean_B_tmp = mean_B
mean_U_tmp = mean_U
mean_B = mean_B_tmp[mean_B_tmp>0]
disp_B = disp_B_tmp[mean_B_tmp>0]
disp_U = disp_U_tmp[mean_U_tmp>0]
mean_U = mean_U_tmp[mean_U_tmp>0]

library(DESeq)
disp_B <- adjustScvForBias( disp_B, 2)
disp_U <- adjustScvForBias( disp_U, 2)


############################################################################################

geneName_B = paste("B",geneName,sep="")
geneName_U = paste("U",geneName,sep="")
data_B = data_div_norm[,1:(ncol(data)/2)]
data_U = data_div_norm[,(ncol(data)/2+1):ncol(data)]
rownames(data_B) = geneName_B
rownames(data_U) = geneName_U
data_BU = rbind(data_B[mean_B_tmp>0,],data_U[mean_U_tmp>0,])

############################################################################################
# Calculate distance matrix
############################################################################################
mat = data_BU
sortLength = nrow(data_BU)

system.time({
smat <- apply(mat, 1, crossprod)
mat1 <- matrix(smat, nrow=sortLength,ncol=sortLength)
mat3 <- tcrossprod(mat)
mat4 <- mat1 + t(mat1) - 2*mat3
 diag(mat4) <- 0
 mat5 <- sqrt(mat4)
 })
comat = mat5

mean_disp_B = cbind(mean_B,disp_B)
mean_disp_U = cbind(mean_U,disp_U)

mean_disp_BU = rbind(mean_disp_B,mean_disp_U)

mean_disp_gene_BU = cbind(mean_disp_BU,as.matrix(rownames(data_BU)))
mean_disp_BU[,1] = log(mean_disp_BU[,1])/100
rownames(mean_disp_BU) = mean_disp_gene_BU[,3]
rm(mat,mat1,mat3,mat4,mat5,sortLength)
data_BU[which(data_BU<.1)]=.00001
data_BU_log = log(data_BU)
dst_BU = density(data_BU_log,n=10000)
############################################################################################
# Fit dispersion
############################################################################################
system.time({
fitted_disp = apply(comat,1,
function(tmp)
{

	tmp_sorted=sort(tmp)
	tmp_names=names(tmp_sorted[1])
	log_tmp = log(mean(data_BU[tmp_names,]))
	log_idx = which(dst_BU$x<log_tmp)
	knn = dst_BU$y[log_idx[length(log_idx)]+1]*nrow(data_BU)/4
	if(knn<100)
		knn = 100
	if(knn>2000)
		knn = 2000
	#tmp_sorted = tmp_sorted[which(tmp_sorted<mean(data_BU[tmp_names,])*0.5)]
	tmp_names = names(tmp_sorted)
	#knn = length(tmp_sorted)
	#knn = 2000
	sortname_tmp = as.matrix(tmp_names)
	mean_disp_gene_i = mean_disp_BU[sortname_tmp,]
	tmp_weight = rep(0,knn)
	tmp_weight_plus_Y = rep(0,knn)
	tmp_weight[1] = 0
	if(knn>1)
	{
		for( i in 2:knn )
		{
		    tmp_dist = sqrt((mean_disp_gene_i[i,1]- mean_disp_gene_i[1,1])^2+ (mean_disp_gene_i[i,2]- mean_disp_gene_i[1,2])^2)
		    if(is.na(tmp_dist) || tmp_dist == 0)
		      tmp_weight[i] = 1
		    else
		      tmp_weight[i] = (1-(abs(( mean_disp_gene_i[i,2]- mean_disp_gene_i[1,2])/tmp_dist))^3)^3
		    tmp_weight_plus_Y[i] = mean_disp_gene_i[i,2]*tmp_weight[i]
		}

		fitted_disp_j = sum(tmp_weight_plus_Y[2:knn])/sum(tmp_weight[2:knn])
	}
	else
		fitted_disp_j = mean_disp_gene_i[2]
	fitted_disp_j

})
})
fitted_disp = as.matrix(fitted_disp)
fitted_disp = t(fitted_disp)

rm(comat)
gc()
save.image()


############################################################################################
# Differential test
############################################################################################

fitted_disp_B = fitted_disp[1:nrow(mean_disp_B)]
fitted_disp_U = fitted_disp[(nrow(mean_disp_BU)-nrow(mean_disp_U)+1):nrow(mean_disp_BU)]
fitted_disp_B_tmp = matrix(0, nrow=nrow(data),ncol=1)
fitted_disp_U_tmp = matrix(0, nrow=nrow(data),ncol=1)

index_B = which(mean_B_tmp>0)
index_U = which(mean_U_tmp>0)

for( i in 1:length(index_B))
{
	fitted_disp_B_tmp[index_B[i]] = fitted_disp_B[i]
}

for( i in 1:length(index_U))
{
	fitted_disp_U_tmp[index_U[i]] = fitted_disp_U[i]	
}

fitted_disp_B_tmp[which(is.na(fitted_disp_B_tmp))]=0
fitted_disp_U_tmp[which(is.na(fitted_disp_U_tmp))]=0

testPval <- MyexactTest(data_div_norm[,1:(ncol(data)/2)],data_div_norm[,(ncol(data)/2+1):ncol(data_div_norm)],fitted_disp_B_tmp,fitted_disp_U_tmp,big.count=1000*ncol(data))
padj = p.adjust( testPval, method="BH" )


###################################################################################################################################################
