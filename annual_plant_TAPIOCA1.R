###############################################################################
#
# Population dynamics of the SPATIALLY IMPLICIT annual plant model for 2 
# species, with calculation of invasion growth rates and coexistence 
# mechanisms
#
# Here, intrinsic fitness, germination, and competition can all be variable
###############################################################################

#=========================================================================
## Load these libraries
#=========================================================================

library(MASS)

#=========================================================================
# Population dynamics implemented numerically with reciprocal invasions
#=========================================================================
#=========================================================================
#Tunable
#=========================================================================

ns=1000 #Lattice width/height -- This should be even
ngens = 8000 #Number of generations (multiple of 4)

#Species traits: 
mFr=matrix(c(2,2))   #Mean reproduction rates
stdFr = matrix(c(sqrt(2),sqrt(2))) #STD deviations
corFr = matrix(c(1, -0.5, -0.5, 1),2,2) #Correlation
sig_Fr= corFr*matrix(c( stdFr[1]^2, stdFr[1]*stdFr[2], stdFr[1]*stdFr[2], 
		stdFr[2]^2),2,2)#Covariance matrix

mgr = matrix(c(0.1,0.1))#Mean germination rates
stdgr = matrix(c(sqrt(0.1),sqrt(0.1))) #STD deviations
corgr = matrix(c(1, -0.1, -0.1, 1),2,2) #Correlation
sig_gr= corgr*matrix(c( stdgr[1]^2, stdgr[1]*stdgr[2], stdgr[1]*stdgr[2], 
		stdgr[2]^2),2,2)#Covariance matrix

m_alphas=matrix( c(0.4,0.4)) #mean competitive impacts
stdalphas = matrix(c(sqrt(0.1),sqrt(0.1))) #STD deviations
coralphas = matrix(c(1, 0.1, 0.1, 1),2,2) #Correlation
sig_alphas= coralphas*matrix(c( stdalphas[1]^2, stdalphas[1]*stdalphas[2], stdalphas[1]*stdalphas[2], 
		stdalphas[2]^2),2,2)#Covariance matrix

sr=0.9 #Survival

#Invasion events (divide into quarters is the easiest)
invasions =c(1, floor(ngens/4), floor(ngens/2), floor(ngens*3/4) )

#=========================================================================
#Internal
#=========================================================================

Fr.t = abs(mvrnorm(ns,mFr,sig_Fr)) #Reproduction spatial series
#Make this an array, joining species spatial map of repro
Fr = array( c(matrix(Fr.t[,1],ngens,ns,byrow=T), matrix(Fr.t[,2],ngens,ns,byrow=T)), dim=c(ngens,ns,2) )

alphas.t = abs(mvrnorm(ns,m_alphas,sig_alphas)) #Competition spatial series
#Standardize values between 0 and 1 based on maximum, then correct distribution.
alphas.t=(alphas.t-(matrix(apply(alphas.t,2,min),ns,2,byrow=T)))/
	((matrix(apply(alphas.t,2,max),ns,2,byrow=T)) -(matrix(apply(alphas.t,2,min),ns,2,byrow=T)))
#Make this an array, joining species spatial map of competition
alphas = array( c(matrix(alphas.t[,1],ngens,ns,byrow=T), matrix(alphas.t[,2],ngens,ns,byrow=T)), dim=c(ngens,ns,2) )


gr.t=abs(mvrnorm(ns,mgr,sig_gr)) #Take the abs value of correlated normal
#Standardize values between 0 and 1 based on maximum,then correct distribution..
gr.t=(gr.t-(matrix(apply(gr.t,2,min),ns,2,byrow=T)))/
	((matrix(apply(gr.t,2,max),ns,2,byrow=T)) -(matrix(apply(gr.t,2,min),ns,2,byrow=T))) 
#Make this an array, joining species spatial map of germination
gr = array( c(matrix(gr.t[,1],ngens,ns,byrow=T), matrix(gr.t[,2],ngens,ns,byrow=T)), dim=c(ngens,ns,2) )


#Population matrixes
#All variation
nrns1_f=matrix(0,ngens,ns)
nrns2_f=matrix(0,ngens,ns)

#Only variation in Fr
nrns1_Fr=matrix(0,ngens,ns)
nrns2_Fr=matrix(0,ngens,ns)

#Only variation in gr
nrns1_gr=matrix(0,ngens,ns)
nrns2_gr=matrix(0,ngens,ns)

#Variation in gr and Fr
nrns1_Fg=matrix(0,ngens,ns)
nrns2_Fg=matrix(0,ngens,ns)

#No variation
nrns1_no=matrix(0,ngens,ns)
nrns2_no=matrix(0,ngens,ns)

nrns2_f[1,]=0.1 #Seed species 2 and allow it to establish as resident
nrns2_Fr[1,]=0.1 #Seed species 2 and allow it to establish as resident
nrns2_gr[1,]=0.1 #Seed species 2 and allow it to establish as resident
nrns2_Fg[1,]=0.1 #Seed species 2 and allow it to establish as resident
nrns2_no[1,]=0.1 #Seed species 2 and allow it to establish as resident



for (n in 1:(ngens-1)) {

	#Invasion: Seed or set to near-zero at the appropriate time steps:
	if (n==invasions[2]) { nrns1_f[n,]=0.0001; nrns1_Fr[n,]=0.0001  
						nrns1_gr[n,]=0.0001; nrns1_Fg[n,]=0.0001 
						nrns1_no[n,]=0.0001}
	#Switch the roles of resident and invader between species 1 and 2
	if (n==invasions[3]) { nrns1_f[n,]=0.1; nrns1_Fr[n,]=0.1
						   nrns1_gr[n,]=0.1; nrns1_Fg[n,]=0.1
						    nrns1_no[n,]=0.1

							nrns2_f[n,]=0; nrns2_Fr[n,]=0
							nrns2_gr[n,]=0; nrns2_Fg[n,]=0 
							nrns2_no[n,]=0   }
	#Second invasion
	if (n==invasions[4]) { nrns2_f[n,]=0.0001; nrns2_Fr[n,]=0.0001
							nrns2_gr[n,]=0.0001; nrns2_Fg[n,]=0.0001 
							nrns2_no[n,]=0.0001}


	#Spatially implicit annual plant model
	nrns1_f[n+1,] = sr*(1-gr[n,,1])*nrns1_f[n,]+Fr[n,,1]*gr[n,,1]*mean(nrns1_f[n,])/(1+alphas[n,,1]*(mean(nrns1_f[n,])*gr[n,,1]+gr[n,,2]*mean(nrns2_f[n,])))
	nrns2_f[n+1,] = sr*(1-gr[n,,2])*nrns2_f[n,]+Fr[n,,2]*gr[n,,2]*mean(nrns2_f[n,])/(1+alphas[n,,2]*(mean(nrns1_f[n,])*gr[n,,1]+gr[n,,2]*mean(nrns2_f[n,])))

	nrns1_Fr[n+1,] = sr*(1-mean(gr[n,,1]))*nrns1_Fr[n,]+Fr[n,,1]*mean(gr[n,,1])*mean(nrns1_Fr[n,])/(1+mean(alphas[n,,1])*(mean(nrns1_Fr[n,])*mean(gr[n,,1])+mean(gr[n,,2])*mean(nrns2_Fr[n,])))
	nrns2_Fr[n+1,] = sr*(1-mean(gr[n,,2]))*nrns2_Fr[n,]+Fr[n,,2]*mean(gr[n,,2])*mean(nrns2_Fr[n,])/(1+mean(alphas[n,,2])*(mean(nrns1_Fr[n,])*mean(gr[n,,1])+mean(gr[n,,2])*mean(nrns2_Fr[n,])))

	nrns1_gr[n+1,] = sr*(1-gr[n,,1])*nrns1_gr[n,]+mean(Fr[n,,1])*gr[n,,1]*mean(nrns1_gr[n,])/(1+mean(alphas[n,,1])*(mean(nrns1_gr[n,])*gr[n,,1]+gr[n,,2]*mean(nrns2_gr[n,])))
	nrns2_gr[n+1,] = sr*(1-gr[n,,2])*nrns2_gr[n,]+mean(Fr[n,,2])*gr[n,,2]*mean(nrns2_gr[n,])/(1+mean(alphas[n,,2])*(mean(nrns1_gr[n,])*gr[n,,1]+gr[n,,2]*mean(nrns2_gr[n,])))

	nrns1_Fg[n+1,] = sr*(1-gr[n,,1])*nrns1_Fg[n,]+Fr[n,,1]*gr[n,,1]*mean(nrns1_Fg[n,])/(1+mean(alphas[n,,1])*(mean(nrns1_Fg[n,])*gr[n,,1]+gr[n,,2]*mean(nrns2_Fg[n,])))
	nrns2_Fg[n+1,] = sr*(1-gr[n,,2])*nrns2_Fg[n,]+Fr[n,,2]*gr[n,,2]*mean(nrns2_Fg[n,])/(1+mean(alphas[n,,2])*(mean(nrns1_Fg[n,])*gr[n,,1]+gr[n,,2]*mean(nrns2_Fg[n,])))

	nrns1_no[n+1,] = sr*(1-mean(gr[n,,1]))*mean(nrns1_no[n,])+mean(Fr[n,,1])*mean(gr[n,,1])*mean(nrns1_no[n,])/(1+mean(alphas[n,,1])*(mean(nrns1_no[n,])*mean(gr[n,,1])+mean(gr[n,,2])*mean(nrns2_no[n,])))
	nrns2_no[n+1,] = sr*(1-mean(gr[n,,2]))*mean(nrns2_no[n,])+mean(Fr[n,,2])*mean(gr[n,,2])*mean(nrns2_no[n,])/(1+mean(alphas[n,,2])*(mean(nrns1_no[n,])*mean(gr[n,,1])+mean(gr[n,,2])*mean(nrns2_no[n,])))


}

lbls1 = c("No Variation","Yield","Germination","Yield+Germ", "Y+G+Alphas")
par(mfrow =c(3,2))
plot(rowMeans(nrns2_no),t="l",main=lbls1[1])
lines(rowMeans(nrns1_no),col="red") 

plot(rowMeans(nrns2_Fr),t="l",main=lbls1[2])
lines(rowMeans(nrns1_Fr),col="red") 

plot(rowMeans(nrns2_gr),t="l",main=lbls1[3])
lines(rowMeans(nrns1_gr),col="red") 

plot(rowMeans(nrns2_Fg),t="l",main=lbls1[4])
lines(rowMeans(nrns1_Fg),col="red") 

plot(rowMeans(nrns2_f),t="l",main=lbls1[5])
lines(rowMeans(nrns1_f),col="red") 

#=========================================================================
# Calculate invasion growth rates and coexistence mechanisms
#=========================================================================

#Fit a linear model to get invasion
a1=invasions[2]
a2=invasions[2]+100
m1=log(rowMeans(nrns1[a1:a2,]))
xx= a1:a2
m1.lm=lm(m1~xx)

#or take the log(gr1[inv+1]/gr1[inv])
nj=invasions[2]
ni=invasions[2]+1

m1d_no = mean(sr*(1-mean(gr[ni,,1]))*mean(nrns1_no[ni,])+mean(Fr[ni,,1])*mean(gr[ni,,1])*mean(nrns1_no[ni,])/(1+mean(alphas[ni,,1])*(mean(nrns1_no[ni,])*mean(gr[ni,,1])+mean(gr[ni,,2])*mean(nrns2_no[ni,]))))/
mean(sr*(1-mean(gr[nj,,1]))*mean(nrns1_no[nj,])+mean(Fr[nj,,1])*mean(gr[nj,,1])*mean(nrns1_no[nj,])/(1+mean(alphas[nj,,1])*(mean(nrns1_no[nj,])*mean(gr[nj,,1])+mean(gr[nj,,2])*mean(nrns2_no[nj,])))) 

# m1d_no2 = mean(sr*(1-rowMeans(gr[r1,,1]))*rowMeans(nrns1_no[r1,])+rowMeans(Fr[r1,,1])*rowMeans(gr[r1,,1])*rowMeans(nrns1_no[r1,])/(1+rowMeans(alphas[r1,,1])*(rowMeans(nrns1_no[r1,])*rowMeans(gr[r1,,1])+rowMeans(gr[r1,,2])*rowMeans(nrns2_no[r1,]))))/
# mean(sr*(1-rowMeans(gr[(r1-1),,1]))*rowMeans(nrns1_no[(r1-1),])+rowMeans(Fr[(r1-1),,1])*rowMeans(gr[(r1-1),,1])*rowMeans(nrns1_no[(r1-1),])/(1+rowMeans(alphas[(r1-1),,1])*(rowMeans(nrns1_no[(r1-1),])*rowMeans(gr[(r1-1),,1])+rowMeans(gr[(r1-1),,2])*rowMeans(nrns2_no[(r1-1),])))) 


m1d_Fr = mean(sr*(1-mean(gr[ni,,1]))*(nrns1_Fr[ni,])+(Fr[ni,,1])*mean(gr[ni,,1])*mean(nrns1_Fr[ni,])/(1+mean(alphas[ni,,1])*(mean(nrns1_Fr[ni,])*mean(gr[ni,,1])+mean(gr[ni,,2])*mean(nrns2_Fr[ni,]))))/
mean(sr*(1-mean(gr[nj,,1]))*mean(nrns1_Fr[nj,])+(Fr[nj,,1])*mean(gr[nj,,1])*mean(nrns1_Fr[nj,])/(1+mean(alphas[nj,,1])*(mean(nrns1_Fr[nj,])*mean(gr[nj,,1])+mean(gr[nj,,2])*mean(nrns2_Fr[nj,])))) 

m1d_gr = mean (sr*(1-(gr[ni,,1]))*(nrns1_gr[ni,])+mean(Fr[ni,,1])*(gr[ni,,1])*mean(nrns1_gr[ni,])/(1+mean(alphas[ni,,1])*(mean(nrns1_gr[ni,])*(gr[ni,,1])+(gr[ni,,2])*mean(nrns2_gr[ni,]))))/
mean(sr*(1-(gr[nj,,1]))*mean(nrns1_gr[nj,])+mean(Fr[nj,,1])*(gr[nj,,1])*mean(nrns1_gr[nj,])/(1+mean(alphas[nj,,1])*(mean(nrns1_gr[nj,])*(gr[nj,,1])+(gr[nj,,2])*mean(nrns2_gr[nj,])))) 

m1d_Fg = mean(sr*(1-(gr[ni,,1]))*(nrns1_Fg[ni,])+(Fr[ni,,1])*(gr[ni,,1])*mean(nrns1_Fg[ni,])/(1+mean(alphas[ni,,1])*(mean(nrns1_Fg[ni,])*(gr[ni,,1])+(gr[ni,,2])*mean(nrns2_Fg[ni,]))))/
mean(sr*(1-(gr[nj,,1]))*mean(nrns1_Fg[nj,])+(Fr[nj,,1])*(gr[nj,,1])*mean(nrns1_Fg[nj,])/(1+mean(alphas[nj,,1])*(mean(nrns1_Fg[nj,])*(gr[nj,,1])+(gr[nj,,2])*mean(nrns2_Fg[nj,])))) 

m1d_f = mean(sr*(1-(gr[ni,,1]))*(nrns1_f[ni,])+(Fr[ni,,1])*(gr[ni,,1])*mean(nrns1_f[ni,])/(1+(alphas[ni,,1])*(mean(nrns1_f[ni,])*(gr[ni,,1])+(gr[ni,,2])*mean(nrns2_f[ni,]))))/
mean(sr*(1-(gr[nj,,1]))*(nrns1_f[nj,])+(Fr[nj,,1])*(gr[nj,,1])*mean(nrns1_f[nj,])/(1+(alphas[nj,,1])*(mean(nrns1_f[nj,])*(gr[nj,,1])+(gr[nj,,2])*mean(nrns2_f[nj,])))) 

m1d_f = mean(sr*(1-(gr[ni,,1]))*(nrns1_f[ni,])+(Fr[ni,,1])*(gr[ni,,1])*mean(nrns1_f[ni,])/(1+(alphas[ni,,1])*(mean(nrns1_f[ni,])*(gr[ni,,1])+(gr[ni,,2])*mean(nrns2_f[ni,]))))/
mean(sr*(1-(gr[nj,,1]))*(nrns1_f[nj,])+(Fr[nj,,1])*(gr[nj,,1])*mean(nrns1_f[nj,])/(1+(alphas[nj,,1])*(mean(nrns1_f[nj,])*(gr[nj,,1])+(gr[nj,,2])*mean(nrns2_f[nj,])))) 


m1d=c(m1d_no,m1d_Fr,m1d_gr,m1d_Fg,m1d_f)

inv_names = c("None", "Yld", "Germ", "Y+G", "Y+G+A")
plot(factor(inv_names),m1d )

#Approximation 1, without the log(E) and log(C) or log gr
r1 = (nj-500):nj #Region where resident is equilibrated

#Mean terms
gim = mean(gr[nj,,1])
grm = mean(gr[nj,,2])
res_eq = mean(nrns2_f[nj,])
#res_eq = mean(nrns2_no[nj,]) #Why does it work with this value??
sm = mean(alphas[ni,,1])
Yix = mean(Fr[ni,,1])
Do=grm*res_eq*sm+1

#Variance terms
ENS = var(alphas[ni,,1])
ENJ = var(nrns2_f[nj,])
EGJ = var(gr[nj,,2])

#Covariance terms
EYSJ = cov(Fr[ni,,1],alphas[ni,,1])
EYGJ = cov(Fr[ni,,1],gr[ni,,2])
EYGI = cov(Fr[ni,,1],gr[ni,,1])
EGSJ = cov(gr[ni,,2],alphas[ni,,1])
EGSI = cov(gr[ni,,1],alphas[ni,,1])
EGIJ = cov(gr[ni,,1],gr[ni,,2])

#Invader growth rate to check the approximation
gr1=-(gim*grm*res_eq*EYSJ)/Do^2-(gim*res_eq*sm*EYGJ)/Do^2+EYGI/Do+
((2*gim*grm^2*res_eq^2*Yix*ENS)/Do^3+(2*gim*res_eq^2*sm^2*Yix*EGJ)/Do^3)/2+
((2*gim*grm*res_eq^2*sm*Yix)/Do^3-(gim*res_eq*Yix)/Do^2)*EGSJ-(grm*res_eq*Yix*EGSI)/Do^2-(res_eq*sm*Yix*EGIJ)/Do^2+
(gim*Yix)/Do+(1-gim)*sr
