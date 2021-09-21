library(MASS)
library(lsei)
library(NMOF)
library(magic)
library(matrixcalc)
########
ONS<-F
########
if(ONS==T) {
dm<-read.table("C:\\Users\\steve\\Desktop\\Male Deaths.csv",header=T,sep=",")
dm<-as.matrix(dm[,-1])
colnames(dm)<-1961:2017
rownames(dm)<-0:105
df<-read.table("C:\\Users\\steve\\Desktop\\Female Deaths.csv",header=T,sep=",")
df<-as.matrix(df[,-1])
colnames(df)<-1961:2017
rownames(df)<-0:105
Em<-read.table("C:\\Users\\steve\\Desktop\\Male Pop.csv",header=T,sep=",")
Em<-as.matrix(Em[,-1])
colnames(Em)<-1961:2017
rownames(Em)<-0:105
Ef<-read.table("C:\\Users\\steve\\Desktop\\Female Pop.csv",header=T,sep=",")
Ef<-as.matrix(Ef[,-1])
colnames(Ef)<-1961:2017
rownames(Ef)<-0:105
} else {
HMD<-read.table("C:\\Users\\steve\\Desktop\\Deaths_1x1 2016.csv",header=T,sep=",")
dm<-matrix(HMD$Male,111,176)
df<-matrix(HMD$Female,111,176)
Em<-matrix(HMD$MalePop,111,176)
Ef<-matrix(HMD$FemalePop,111,176)
colnames(dm)<-colnames(df)<-colnames(Em)<-colnames(Ef)<-1841:2016
rownames(dm)<-rownames(df)<-rownames(Em)<-rownames(Ef)<-0:110
dm<-dm[,(1961-1841+1):(2016-1841+1)]
df<-df[,(1961-1841+1):(2016-1841+1)]
Em<-Em[,(1961-1841+1):(2016-1841+1)]
Ef<-Ef[,(1961-1841+1):(2016-1841+1)]
}
HMD.Scot<-read.table("C:\\Users\\steve\\Desktop\\Deaths_1x1 2016 SC.csv",header=T,sep=",")
dm.Scot<-matrix(round(HMD.Scot$Male),111,162)
df.Scot<-matrix(round(HMD.Scot$Female),111,162)
Em.Scot<-matrix(HMD.Scot$MalePop,111,162)
Ef.Scot<-matrix(HMD.Scot$FemalePop,111,162)
colnames(dm.Scot)<-colnames(df.Scot)<-colnames(Em.Scot)<-colnames(Ef.Scot)<-1855:2016
rownames(dm.Scot)<-rownames(df.Scot)<-rownames(Em.Scot)<-rownames(Ef.Scot)<-0:110

no.basis=40;extra.knots=0
last.age=104;no.years=56;start.year=1961;end.year=start.year+no.years-1
last.age.Scot=99;no.years.Scot=no.years;start.year.Scot=1961;end.year.Scot=start.year.Scot+no.years.Scot-1
cohort.omit=6

datam=round(c(dm[2:(last.age+1),1:no.years]))
dataf=round(c(df[2:(last.age+1),1:no.years]))
offsetm=log(c(Em[2:(last.age+1),1:no.years]))
offsetf=log(c(Ef[2:(last.age+1),1:no.years]))
datam.Scot=round(c(dm.Scot[2:(last.age.Scot+1),start.year.Scot:end.year.Scot-1855+1]))
dataf.Scot=round(c(df.Scot[2:(last.age.Scot+1),start.year.Scot:end.year.Scot-1855+1]))
offsetm.Scot=log(c(Em.Scot[2:(last.age.Scot+1),start.year.Scot:end.year.Scot-1855+1]))
offsetf.Scot=log(c(Ef.Scot[2:(last.age.Scot+1),start.year.Scot:end.year.Scot-1855+1]))

knots<-seq(0,1,length=no.basis-2)
dk<-knots[2]-knots[1]	
knots<-c(knots[1]-dk*(3:1),knots,knots[no.basis-2]+dk*(1:(3+extra.knots)))
age<-seq(0,1,length=last.age)
age.Scot<-((1:last.age.Scot)-1)/(last.age-1)
cohort.years<-seq(0,1,length=last.age+no.years-1-cohort.omit)

bspline<-function (x,k,i,m=2) {
if (m==-1) {basis<-as.numeric(x<k[i+1] & x>=k[i])} else {
	z0<-(x-k[i])/(k[i+m+1]-k[i])
	z1<-(k[i+m+2]-x)/(k[i+m+2]-k[i+1])
	basis<-z0*bspline(x,k,i,m-1)+z1*bspline(x,k,i+1,m-1) }
basis
}

timeindex<-((1:no.years)-1)/(no.years-1)
timeindex<-timeindex-mean(timeindex)
timeindex.Scot<-timeindex[(start.year.Scot-start.year+1):(end.year.Scot-start.year+1)]

A<-c()
for(j in 1:no.basis) {
A<-cbind(A,bspline(age,knots,j))
}

A.Scot<-c()
for(j in 1:no.basis) {
A.Scot<-cbind(A.Scot,bspline(age.Scot,knots,j))
}

DX<-rep(1,no.years)%x%A
DX.Scot<-rep(1,no.years.Scot)%x%A.Scot
DX.period<-timeindex%x%A
DX.period.Scot<-timeindex.Scot%x%A.Scot

kappa.trans<-rbind(rep(1,no.years),timeindex,cbind(matrix(0,no.years-2,2),diag(no.years-2)))
inv.kappa.trans<-solve(kappa.trans)
kappa.mat<-rbind(t(inv.kappa.trans[1,3:no.years])%x%rep(1,last.age),t(inv.kappa.trans[2,3:no.years])%x%rep(1,last.age),diag(no.years-2)%x%rep(1,last.age))
kappa.trans.Scot<-rbind(rep(1,no.years.Scot),timeindex.Scot,cbind(matrix(0,no.years.Scot-2,2),diag(no.years.Scot-2)))
inv.kappa.trans.Scot<-solve(kappa.trans.Scot)
kappa.mat.Scot<-rbind(t(solve(kappa.trans.Scot)[1,3:no.years.Scot])%x%rep(1,last.age.Scot),t(solve(kappa.trans.Scot)[2,3:no.years.Scot])%x%rep(1,last.age.Scot),diag(no.years.Scot-2)%x%rep(1,last.age.Scot))

A2<-c()
for(j in 1:no.basis) {
A2<-cbind(A2,bspline(cohort.years,knots,j))
}

CX.constraints<-rbind(A2[1,],apply(A2,2,sum),A2[nrow(A2),])
C.Q<-qr.Q(qr(t(CX.constraints)),complete=T)
C.Q.trans<-C.Q[,-(1:3)]
A2.reduced<-A2%*%C.Q.trans
A2.reduced<-rbind(matrix(0,cohort.omit,no.basis-3),A2.reduced)

A2.Scot<-A2[((start.year.Scot-last.age.Scot)-(start.year-last.age)+1):((start.year.Scot-last.age.Scot)-(start.year-last.age)+no.years.Scot+last.age.Scot-1-cohort.omit),]
A2.Scot<-A2.Scot[,colSums(A2.Scot)!=0]
CX.constraints.Scot<-rbind(A2.Scot[1,],apply(A2.Scot,2,sum),A2.Scot[nrow(A2.Scot),])
C.Q.Scot<-qr.Q(qr(t(CX.constraints.Scot)),complete=T)
C.Q.trans.Scot<-C.Q.Scot[,-(1:3)]
A2.reduced.Scot<-A2.Scot%*%C.Q.trans.Scot
A2.reduced.Scot<-rbind(matrix(0,cohort.omit,ncol(A2.Scot)-3),A2.reduced.Scot)

CX<-c()
for(i in 1:no.years){
CX<-rbind(CX,A2.reduced[(i-1)+(last.age:1),])
}
CX2<-c()
for(i in 1:no.years.Scot){
CX2<-rbind(CX2,A2.reduced.Scot[(i-1)+(last.age.Scot:1),])
}

DX.cohort<-CX
DX.cohort.Scot<-CX2

omit.index<-c();for(i in 1:cohort.omit){omit.index<-c(omit.index,last.age*(i-1)+(last.age-cohort.omit+i):last.age)}
omit.index.Scot<-c();for(i in 1:cohort.omit){omit.index.Scot<-c(omit.index.Scot,last.age.Scot*(i-1)+(last.age.Scot-cohort.omit+i):last.age.Scot)}

DX.all<-cbind(adiag(DX[-omit.index,],DX[-omit.index,],DX.Scot[-omit.index.Scot,],DX.Scot[-omit.index.Scot,]),adiag(DX.period[-omit.index,],DX.period[-omit.index,],DX.period.Scot[-omit.index.Scot,],DX.period.Scot[-omit.index.Scot,]),adiag(kappa.mat[-omit.index,],kappa.mat[-omit.index,],kappa.mat.Scot[-omit.index.Scot,],kappa.mat.Scot[-omit.index.Scot,]),adiag(DX.cohort[-omit.index,],DX.cohort[-omit.index,],DX.cohort.Scot[-omit.index.Scot,],DX.cohort.Scot[-omit.index.Scot,]))
DX.all.Eng<-cbind(adiag(DX[-omit.index,],DX[-omit.index,]),adiag(DX.period[-omit.index,],DX.period[-omit.index,]),adiag(kappa.mat[-omit.index,],kappa.mat[-omit.index,]),adiag(DX.cohort[-omit.index,],DX.cohort[-omit.index,]))
DX.all.Scot<-cbind(adiag(DX.Scot[-omit.index.Scot,],DX.Scot[-omit.index.Scot,]),adiag(DX.period.Scot[-omit.index.Scot,],DX.period.Scot[-omit.index.Scot,]),adiag(kappa.mat.Scot[-omit.index.Scot,],kappa.mat.Scot[-omit.index.Scot,]),adiag(DX.cohort.Scot[-omit.index.Scot,],DX.cohort.Scot[-omit.index.Scot,]))
DX.sex<-cbind(adiag(DX[-omit.index,],DX.Scot[-omit.index.Scot,]),adiag(DX.period[-omit.index,],DX.period.Scot[-omit.index.Scot,]),adiag(kappa.mat[-omit.index,],kappa.mat.Scot[-omit.index.Scot,]),adiag(DX.cohort[-omit.index,],DX.cohort.Scot[-omit.index.Scot,]))
DX.all.single<-cbind(DX,DX.period,kappa.mat,DX.cohort)[-omit.index,]
DX.all.single.Scot<-cbind(DX.Scot,DX.period.Scot,kappa.mat.Scot,DX.cohort.Scot)[-omit.index.Scot,]

DX.all.single.nc<-cbind(DX,DX.period,kappa.mat)[-omit.index,]
DX.all.single.Scot.nc<-cbind(DX.Scot,DX.period.Scot,kappa.mat.Scot)[-omit.index.Scot,]
DX.all.Eng.nc<-cbind(adiag(DX[-omit.index,],DX[-omit.index,]),adiag(DX.period[-omit.index,],DX.period[-omit.index,]),adiag(kappa.mat[-omit.index,],kappa.mat[-omit.index,]))
DX.all.Scot.nc<-cbind(adiag(DX.Scot[-omit.index.Scot,],DX.Scot[-omit.index.Scot,]),adiag(DX.period.Scot[-omit.index.Scot,],DX.period.Scot[-omit.index.Scot,]),adiag(kappa.mat.Scot[-omit.index.Scot,],kappa.mat.Scot[-omit.index.Scot,]))
DX.all.Eng.nc<-cbind(adiag(DX[-omit.index,],DX[-omit.index,]),adiag(DX.period[-omit.index,],DX.period[-omit.index,]),adiag(kappa.mat[-omit.index,],kappa.mat[-omit.index,]))
DX.all.Scot.nc<-cbind(adiag(DX.Scot[-omit.index.Scot,],DX.Scot[-omit.index.Scot,]),adiag(DX.period.Scot[-omit.index.Scot,],DX.period.Scot[-omit.index.Scot,]),adiag(kappa.mat.Scot[-omit.index.Scot,],kappa.mat.Scot[-omit.index.Scot,]))

P<-diff(diag(no.basis+extra.knots),differences=2)
PD<-t(c(1,-1))%x%diag(no.basis+extra.knots)
PD1<-t(c(1,-1))%x%diff(diag(no.basis+extra.knots),differences=1)
PD1.joint<-t(c(1,0,-1))%x%diff(diag(no.basis+extra.knots),differences=1)
#PC<-(diff(diag(no.basis),differences=2)%*%solve(CX.constraints))[,-c(1,2,no.basis)]
PC<-diff(diag(no.basis),differences=2)%*%C.Q.trans
PC.Scot<-diff(diag(ncol(A2.Scot)),differences=2)%*%C.Q.trans.Scot

P.index.single<-list(1,
			no.basis+1,
			no.basis*2+no.years-2+1)

P.index.PD1<-list(c(1,no.basis+1),
			c(1),
			c(no.basis*2+1,no.basis*3+1),
			c(no.basis*2+1),
			c(no.basis*4+(no.years-2)*2+1,no.basis*4+(no.years-2)*2+no.basis-3+1))
P.index.PD.PD1<-list(c(1,no.basis+1),
			c(1),
			c(1),
			c(no.basis*2+1,no.basis*3+1),
			c(no.basis*2+1),
			c(no.basis*2+1),
			c(no.basis*4+(no.years-2)*2+1,no.basis*4+(no.years-2)*2+no.basis-3+1))
P.index.sex<-list(c(1,no.basis+1),
			c(no.basis*2+1,no.basis*3+1),
			c(no.basis*4+(no.years-2)*2+1,no.basis*4+(no.years-2)*2+no.basis-3+1))
P.index.Eng<-list(1,
			no.basis+1,
			1,
			no.basis*2+1,
			no.basis*3+1,
			no.basis*2+1,
			no.basis*4+(no.years-2)*2+1,
			no.basis*4+(no.years-2)*2+no.basis-3+1)
P.index.Scot<-list(1,
			no.basis+1,
			1,
			no.basis*2+1,
			no.basis*3+1,
			no.basis*2+1,
			no.basis*4+(no.years-2)*2+1,
			no.basis*4+(no.years-2)*2+ncol(A2.reduced.Scot)+1)
P.index.Eng.nc<-list(1,
			no.basis+1,
			1,
			no.basis*2+1,
			no.basis*3+1,
			no.basis*2+1)
P.index.Scot.nc<-list(1,
			no.basis+1,
			1,
			no.basis*2+1,
			no.basis*3+1,
			no.basis*2+1)
P.index.Eng.Scot<-list(c(1,(no.basis+extra.knots)*2+1),
			c((no.basis+extra.knots)+1,(no.basis+extra.knots)*3+1),
			c(1,(no.basis+extra.knots)*2+1),
			c(1,(no.basis+extra.knots)+1),
			c((no.basis+extra.knots)*4+1,(no.basis+extra.knots)*6+1),
			c((no.basis+extra.knots)*5+1,(no.basis+extra.knots)*7+1),
			c((no.basis+extra.knots)*4+1,(no.basis+extra.knots)*6+1),
			c((no.basis+extra.knots)*4+1,(no.basis+extra.knots)*5+1),
			c((no.basis+extra.knots)*8+(no.years-2)*4+1,(no.basis+extra.knots)*8+(no.years-2)*4+(no.basis-3)*2+1),
			c((no.basis+extra.knots)*8+(no.years-2)*4+no.basis-3+1,(no.basis+extra.knots)*8+(no.years-2)*4+(no.basis-3)*2+ncol(A2.reduced.Scot)+1))
P.list.single<-list(list(P),list(P),list(PC))
P.list.single.Scot<-list(list(P),list(P),list(PC.Scot))
P.list.Eng<-list(list(P),list(P),list(PD[-(1:8),]),list(P),list(P),list(PD[-(1:8),]),list(PC),list(PC))
P.list.Scot<-list(list(P),list(P),list(PD[-(1:8),]),list(P),list(P),list(PD[-(1:8),]),list(PC.Scot),list(PC.Scot))
P.list.PD1=list(list(P,P),list(PD1),list(P,P),list(PD1),list(PC,PC.Scot))
P.list.PD.PD1=list(list(P,P),list(PD),list(PD1),list(P,P),list(PD),list(PD1),list(PC,PC.Scot))
P.list.sex=list(list(P,P),list(P,P),list(PC,PC.Scot))
P.list.Eng.nc<-list(list(P),list(P),list(PD[-(1:8),]),list(P),list(P),list(PD[-(1:8),]))
P.list.Scot.nc<-list(list(P),list(P),list(PD[-(1:8),]),list(P),list(P),list(PD[-(1:8),]))

build.penalty<-function(P,index,total.col){
B<-c()
for(i in 1:length(P)) {
B<-rbind(B,cbind(matrix(0,nrow(P[[i]]),index[i]-1),P[[i]],matrix(0,nrow(P[[i]]),total.col-ncol(P[[i]])-index[i]+1)))
	}
B
}

bic.function<-function(lambda,X,y,P,exprate.index,H=NULL) {
cat(lambda,"\n")

X.QR<-qr(X)
QR.pivot<-X.QR$pivot
Q<-qr.Q(X.QR)
R<-qr.R(X.QR)
total.col<-ncol(X)
P.weighted<-list()
for(i in 1:length(P)){
j<-sum(exprate.index[0:(i-1)])+(i-1)+1
P[[i]]<-lapply(P[[i]],function(X){X<-X[,QR.pivot]})
if(exprate.index[i]==1) {lambda.weights<-exp(lambda[j]+lambda[j+1]*seq(0,1,length=max(sapply(P[[i]],nrow))));P.weighted[[i]]<-lapply(P[[i]],function(X){as.numeric(lambda.weights[1:nrow(X)]^0.5)*X})} else {lambda.weights<-exp(lambda[j]);P.weighted[[i]]<-lapply(P[[i]],function(X){as.numeric(lambda.weights^0.5)*X})}
}

B<-c()
for(i in 1:length(P.weighted)){
	for(j in 1:length(P.weighted[[i]])){
B<-rbind(B,P.weighted[[i]][[j]])
	}
}


SVD<-svd(rbind(R,B,H[,QR.pivot]))
#####################lowest svd d
#svd.index<-which(SVD$d > max(SVD$d)*sqrt(.Machine$double.eps))
#####################
svd.index<-which(SVD$d > max(SVD$d)*1e-7)
V<-SVD$v[,svd.index]
U1<-SVD$u[1:ncol(X),svd.index]
D<-SVD$d[svd.index]
M<-rep(list(0),sum(exprate.index)+length(exprate.index))
for(i in 1:length(P)){
j<-sum(exprate.index[0:(i-1)])+(i-1)+1
for(k in 1:length(P.weighted[[i]])){
temp.index<-apply(P.weighted[[i]][[k]],2,function(x){!all(x==0)})
M[[j]]<-M[[j]]+crossprod(P.weighted[[i]][[k]][,temp.index]%*%V[temp.index,]%*%diag(1/D))
if (exprate.index[i]==1) {
M[[j+1]]<-M[[j+1]]+crossprod(diag(sqrt(seq(0,1,length=max(sapply(P[[i]],nrow)))[1:nrow(P[[i]][[k]])]))%*%P.weighted[[i]][[k]][,temp.index]%*%V[temp.index,]%*%diag(1/D))
		}
	}
}
K<-lapply(M,function(X){X%*%t(U1)%*%U1})

err<-as.vector(y-(Q%*%(U1%*%(t(U1)%*%(t(Q)%*%y)))))
y.star<-t(U1)%*%(t(Q)%*%y)
err.star<-t(U1)%*%(t(Q)%*%err)
sigma.2<-as.numeric(crossprod(err)/length(y))
dsig.dlambda<-trA.mat<-c()
for(i in 1:length(lambda)){
dsig.dlambda[i]<-2/length(y)*(t(err.star)%*%(M[[i]]%*%y.star))
trA.mat[i]<-sum(diag(K[[i]]))
}
d2sig.dlambda.base<-tr2A.mat.base<-matrix(0,nrow=length(lambda),ncol=length(lambda))
for(i in 1:length(lambda)) {
	for(j in i:length(lambda)) {
d2sig.dlambda.base[i,j]<--2/length(y)*(t(err.star)%*%((M[[i]]%*%M[[j]]+M[[j]]%*%M[[i]])%*%y.star)-t(t(K[[i]])%*%y.star)%*%(M[[j]]%*%y.star))
d2sig.dlambda.base[j,i]<-d2sig.dlambda.base[i,j]
tr2A.mat.base[i,j]<-sum(diag(M[[i]]%*%K[[j]]))+sum(diag(M[[j]]%*%K[[i]]))
tr2A.mat.base[j,i]<-tr2A.mat.base[i,j]
	}
}
d2sig.dlambda.extra<-tr2A.mat.extra<-matrix(0,nrow=length(lambda),ncol=length(lambda))
for(i in 1:length(exprate.index)){
j<-sum(exprate.index[0:(i-1)])+(i-1)+1
d2sig.dlambda.extra[j,j]<-dsig.dlambda[j]
tr2A.mat.extra[j,j]<-trA.mat[j]
if(exprate.index[i]==1) {
for(k in 1:length(P.weighted[[i]])){
temp.index<-apply(P.weighted[[i]][[k]],2,function(x){!all(x==0)})
d2sig.dlambda.extra[j+1,j+1]<-d2sig.dlambda.extra[j+1,j+1]+2/length(y)*t(err.star)%*%(crossprod(diag(seq(0,1,length=max(sapply(P[[i]],nrow)))[1:nrow(P[[i]][[k]])])%*%P.weighted[[i]][[k]][,temp.index]%*%V[temp.index,]%*%diag(1/D))%*%y.star)
tr2A.mat.extra[j+1,j+1]<-tr2A.mat.extra[j+1,j+1]+sum(diag(crossprod(diag(seq(0,1,length=max(sapply(P[[i]],nrow)))[1:nrow(P[[i]][[k]])])%*%P.weighted[[i]][[k]][,temp.index]%*%V[temp.index,]%*%diag(1/D))%*%t(U1)%*%U1))
		}
d2sig.dlambda.extra[j,j+1]<-dsig.dlambda[j+1]
d2sig.dlambda.extra[j+1,j]<-d2sig.dlambda.extra[j,j+1]
tr2A.mat.extra[j,j+1]<-trA.mat[j+1]
tr2A.mat.extra[j+1,j]<-tr2A.mat.extra[j,j+1]
	}
}

d2sig.dlambda<-d2sig.dlambda.base+d2sig.dlambda.extra
tr2A.mat<-tr2A.mat.base-tr2A.mat.extra
#bic<-length(y)*log(sigma.2)+log(length(y))*sum(diag(tcrossprod(U1)))
#attr(bic,"gradient")<- (length(y)/sigma.2)*dsig.dlambda-log(length(y))*trA.mat
#attr(bic,"hessian")<-(length(y)/sigma.2)*(d2sig.dlambda-(1/sigma.2)*tcrossprod(dsig.dlambda))+log(length(y))*tr2A.mat

bic<-(length(y)*sigma.2+log(length(y))*sum(diag(tcrossprod(U1))))
attr(bic,"gradient")<-(length(y)*dsig.dlambda-log(length(y))*trA.mat)
attr(bic,"hessian")<-(length(y)*d2sig.dlambda+log(length(y))*tr2A.mat)
bic

#alpha<-length(y)*sigma.2
#delta<-length(y)-sum(diag(tcrossprod(U1)))
#gcv<-length(y)*(length(y)*alpha/delta^2)
#attr(gcv,"gradient")<-length(y)*(length(y)^2/delta^2*dsig.dlambda-2*length(y)*alpha/delta^3*trA.mat)
#attr(gcv,"hessian")<-length(y)*(-2*length(y)^2/delta^3*(dsig.dlambda%*%t(trA.mat)+trA.mat%*%t(dsig.dlambda))+length(y)/delta^2*d2sig.dlambda+6*length(y)*alpha/delta^4*(trA.mat%*%t(trA.mat))-2*length(y)*alpha/delta^3*d2sig.dlambda)
#gcv
}


newton.raphson.bic<-function(data,mu=data,offset,DXX,P.list,P.index,exp.index,lambda,H=NULL,fnval,stepmax,steptol,gradtol){
mu<-ifelse(mu==0,1e-5,mu)
dev<-ll.sat<-sum(dpois(data,data,log=TRUE))
eta<-log(mu)
for(i in 1:length(P.list)){
	for(j in 1:length(P.list[[i]])){
P.list[[i]][[j]]<-build.penalty(P.list[[i]][j],index=P.index[[i]][j],total.col=ncol(DXX))
}
}
converged.magic<-FALSE
while(!converged.magic) {
	z<-(data-mu)/mu+eta-offset
	w<-mu
	min.bic<-nlm(f=bic.function,p=lambda,X=as.numeric(w^0.5)*DXX,y=as.numeric(w^0.5)*z,P=P.list,exprate.index=exp.index,H=H,print.level=2,gradtol=gradtol,check.analyticals=F,iterlim=1000,typsize=abs(lambda),stepmax=stepmax,fscale=abs(fnval),steptol=steptol)
if(min.bic$code!=1 & min.bic$code!=2 & min.bic$code!=3) stop("nlm code not 1 or 2 or 3")	
if(min.bic$code==3) cat("nlm code 3")
	lambda<-min.bic$estimate
	fnval<-min.bic$minimum
P.weighted<-list()
for(i in 1:length(P.list)){
j<-sum(exp.index[0:(i-1)])+(i-1)+1
if(exp.index[i]==1) {lambda.weights<-exp(lambda[j]+lambda[j+1]*seq(0,1,length=max(sapply(P.list[[i]],nrow))));P.weighted[[i]]<-lapply(P.list[[i]],function(X){as.numeric(lambda.weights[1:nrow(X)]^0.5)*X})} else {lambda.weights<-exp(lambda[j]);P.weighted[[i]]<-lapply(P.list[[i]],function(X){as.numeric(lambda.weights^0.5)*X})}
}
B<-c()
for(i in 1:length(P.weighted)){
	for(j in 1:length(P.weighted[[i]])){
B<-rbind(B,P.weighted[[i]][[j]])
	}
}
fit<-lm(c(as.numeric(w^0.5)*z,rep(0,nrow(rbind(H,B))))~rbind(as.numeric(w^0.5)*DXX,H,B)-1)
	eta<-DXX%*%coef(fit)+offset
	mu<-exp(eta)
	old.dev<-dev
	dev<-2*(ll.sat-sum(dpois(data,mu,log=TRUE)))
	if(abs(dev-old.dev)<1e-6*dev) converged.magic<-TRUE
}

list(b=coef(fit),sp=lambda,bic=min.bic$minimum,gradient=min.bic$gradient,hessian=min.bic$hessian,B=B,fv=mu)
}


#####OPTIM
lambda.optim<-function(lambda,y,XX,P.list,P.index,exp.index,H=NULL){
P.weighted<-list()
for(i in 1:length(P.list)){
j<-sum(exp.index[0:(i-1)])+(i-1)+1
P.weighted[[i]]<-lapply(P.list[[i]],function(X){if(exp.index[i]==1) lambda.weights<-exp(lambda[j]+lambda[j+1]*seq(0,1,length=nrow(X))) else lambda.weights<-exp(lambda[j]); as.numeric(lambda.weights^0.5)*X})
}
B.list<-list()
for(i in 1:length(P.weighted)){
B.list<-c(B.list,P.weighted[[i]])
}

B<-build.penalty(B.list,unlist(P.index),ncol(XX))

fit<-lm(c(y,rep(0,nrow(rbind(H,B))))~rbind(XX,H,B)-1)
length(y)*log(crossprod(y-fitted(fit)[1:length(y)])/length(y))+log(length(y))*sum(influence(fit)$hat[1:length(y)])
}

newton.raphson.bic.grid<-function(lambda,data,offset,DXX,P.list,exp.index,H=NULL,typsize,fscale,stepmax){
mu<-ifelse(data==0,1e-5,data)
dev<-ll.sat<-sum(dpois(data,data,log=TRUE))
eta<-log(mu)
exprate.index<-exp.index
converged.magic<-FALSE
while(!converged.magic) {
	z<-(data-mu)/mu+eta-offset
	w<-mu
P.weighted<-list()
for(i in 1:length(P.list)){
j<-sum(exprate.index[0:(i-1)])+(i-1)+1
P.weighted[[i]]<-lapply(P.list[[i]],function(X){if(exprate.index[i]==1) lambda.weights<-exp(lambda[j]+lambda[j+1]*seq(0,1,length=nrow(X))) else lambda.weights<-exp(lambda[j]); as.numeric(lambda.weights^0.5)*X})
}

B<-c()
for(i in 1:length(P.weighted)){
	for(j in 1:length(P.weighted[[i]])){
B<-rbind(B,P.weighted[[i]][[j]])
	}
}

fit<-lm(c(as.numeric(w^0.5)*z,rep(0,nrow(rbind(H,B))))~rbind(as.numeric(w^0.5)*DXX,H,B)-1)
	eta<-DXX%*%coef(fit)+offset
	mu<-exp(eta)
	old.dev<-devs
	dev<-2*(ll.sat-sum(dpois(data,mu,log=TRUE)))
	if(abs(dev-old.dev)<1e-6*dev) converged.magic<-TRUE
}
length(z)*log(crossprod(as.numeric(w^0.5)*z-fitted(fit)[1:length(z)])/length(z))+log(length(z))*sum(influence(fit)$hat[1:length(z)])
}

Eng<-newton.raphson.bic(data=datam[-omit.index],offset=offsetm[-omit.index],DXX=DX.all.single,P.list=P.list.single,P.index=P.index.single,exp.index=c(1,0,0),lambda=c(2,4,6,6),fnval=1e4,stepmax=10,steptol=1e-12,gradtol=1e-8)
Scot<-newton.raphson.bic(data=datam.Scot[-omit.index.Scot],offset=offsetm.Scot[-omit.index.Scot],DXX=DX.all.single.Scot,P.list=P.list.single.Scot,P.index=P.index.single,exp.index=c(1,0,0),lambda=c(2,3,6,6),fnval=1e4,stepmax=10,steptol=1e-12,gradtol=1e-8)
Engf<-newton.raphson.bic(data=dataf[-omit.index],offset=offsetf[-omit.index],DXX=DX.all.single,P.list=P.list.single,P.index=P.index.single,exp.index=c(1,0,0),lambda=c(2,4,6,6),fnval=1e4,stepmax=10,steptol=1e-12,gradtol=1e-8)
Scotf<-newton.raphson.bic(data=dataf.Scot[-omit.index.Scot],offset=offsetf.Scot[-omit.index.Scot],DXX=DX.all.single.Scot,P.list=P.list.single.Scot,P.index=P.index.single,exp.index=c(1,0,0),lambda=c(2,4,6,6),fnval=1e4,stepmax=10,steptol=1e-12,gradtol=1e-8)

Eng.46<-newton.raphson.bic(data=datam[-omit.index],offset=offsetm[-omit.index],DXX=DX.all.single,P.list=P.list.single,P.index=P.index.single,exp.index=c(1,0,0),lambda=c(2,4,6,6),fnval=1e4,stepmax=10,steptol=1e-12,gradtol=1e-8)
Scot.46<-newton.raphson.bic(data=datam.Scot[-omit.index.Scot],offset=offsetm.Scot[-omit.index.Scot],DXX=DX.all.single.Scot,P.list=P.list.single.Scot,P.index=P.index.single,exp.index=c(1,0,0),lambda=c(2,3,6,6),fnval=1e4,stepmax=10,steptol=1e-12,gradtol=1e-8)
Engf.46<-newton.raphson.bic(data=dataf[-omit.index],offset=offsetf[-omit.index],DXX=DX.all.single,P.list=P.list.single,P.index=P.index.single,exp.index=c(1,0,0),lambda=c(2,4,6,6),fnval=1e4,stepmax=10,steptol=1e-12,gradtol=1e-8)
Scotf.46<-newton.raphson.bic(data=dataf.Scot[-omit.index.Scot],offset=offsetf.Scot[-omit.index.Scot],DXX=DX.all.single.Scot,P.list=P.list.single.Scot,P.index=P.index.single,exp.index=c(1,0,0),lambda=c(2,4,6,6),fnval=1e4,stepmax=10,steptol=1e-12,gradtol=1e-8)

Eng.56<-newton.raphson.bic(data=datam[-omit.index],offset=offsetm[-omit.index],DXX=DX.all.single,P.list=P.list.single,P.index=P.index.single,exp.index=c(1,0,0),lambda=c(2,4,6,6),typsize=c(1,1,1.5,1.5),fscale=1e3,stepmax=10)
Scot.56<-newton.raphson.bic(data=datam.Scot[-omit.index.Scot],offset=offsetm.Scot[-omit.index.Scot],DXX=DX.all.single.Scot,P.list=P.list.single.Scot,P.index=P.index.single,exp.index=c(1,0,0),lambda=c(-1,3,6,6),typsize=c(1,1e-1,1.5,1.5),fscale=1,stepmax=30)
Engf.56<-newton.raphson.bic(data=dataf[-omit.index],offset=offsetf[-omit.index],DXX=DX.all.single,P.list=P.list.single,P.index=P.index.single,exp.index=c(1,0,0),lambda=c(2,4,6,6),typsize=c(1,1,1.5,1.5),fscale=1e3,stepmax=10)
Scotf.56<-newton.raphson.bic(data=dataf.Scot[-omit.index.Scot],offset=offsetf.Scot[-omit.index.Scot],DXX=DX.all.single.Scot,P.list=P.list.single.Scot,P.index=P.index.single,exp.index=c(1,0,0),lambda=c(2,4,6,6),typsize=c(1,1,1.5,1.5),fscale=1e3,stepmax=10)

Eng.nc<-newton.raphson.bic(data=datam[-omit.index],offset=offsetm[-omit.index],DXX=DX.all.single.nc,P.list=P.list.single[-3],P.index=P.index.single[-3],exp.index=c(1,0),lambda=c(2,4,6),typsize=c(1,1,1.5),fscale=1e3,stepmax=10)
Scot.nc<-newton.raphson.bic(data=datam.Scot[-omit.index.Scot],offset=offsetm.Scot[-omit.index.Scot],DXX=DX.all.single.Scot.nc,P.list=P.list.single.Scot[-3],P.index=P.index.single[-3],exp.index=c(1,0),lambda=c(-1,3,6),typsize=c(1,1,1.5),fscale=1,stepmax=30)
Engf.nc<-newton.raphson.bic(data=dataf[-omit.index],offset=offsetf[-omit.index],DXX=DX.all.single.nc,P.list=P.list.single[-3],P.index=P.index.single[-3],exp.index=c(1,0),lambda=c(2,4,6),typsize=c(1,1,1.5),fscale=1e3,stepmax=10)
Scotf.nc<-newton.raphson.bic(data=dataf.Scot[-omit.index.Scot],offset=offsetf.Scot[-omit.index.Scot],DXX=DX.all.single.Scot.nc,P.list=P.list.single.Scot[-3],P.index=P.index.single[-3],exp.index=c(1,0),lambda=c(2,4,6),typsize=c(1,1,1.5),fscale=1e3,stepmax=10)

lambda.weightsm<-exp(Eng$sp[1]+Eng$sp[2]*seq(0,1,length=nrow(P)))
lambda.weightsf<-exp(Engf$sp[1]+Engf$sp[2]*seq(0,1,length=nrow(P)))
lambda.weightsm.Scot<-exp(Scot$sp[1]+Scot$sp[2]*seq(0,1,length=nrow(P)))
lambda.weightsf.Scot<-exp(Scotf$sp[1]+Scotf$sp[2]*seq(0,1,length=nrow(P)))
P.alpha.m<-as.numeric(lambda.weightsm^0.5)*P
P.alpha.f<-as.numeric(lambda.weightsf^0.5)*P
P.alpha.m.Scot<-as.numeric(lambda.weightsm.Scot^0.5)*P
P.alpha.f.Scot<-as.numeric(lambda.weightsf.Scot^0.5)*P
P.beta.m<-as.numeric(exp(Eng$sp[3])^0.5)*P
P.beta.f<-as.numeric(exp(Engf$sp[3])^0.5)*P
P.beta.m.Scot<-as.numeric(exp(Scot$sp[3])^0.5)*P
P.beta.f.Scot<-as.numeric(exp(Scotf$sp[3])^0.5)*P
P.gamma.m<-as.numeric(exp(Eng$sp[4])^0.5)*PC
P.gamma.f<-as.numeric(exp(Engf$sp[4])^0.5)*PC
P.gamma.m.Scot<-as.numeric(exp(Scot$sp[4])^0.5)*PC.Scot
P.gamma.f.Scot<-as.numeric(exp(Scotf$sp[4])^0.5)*PC.Scot
P.offset.list.Eng<-list(P.alpha.m,P.alpha.f,P.beta.m,P.beta.f,P.gamma.m,P.gamma.f)
P.offset.list.Scot<-list(P.alpha.m.Scot,P.alpha.f.Scot,P.beta.m.Scot,P.beta.f.Scot,P.gamma.m.Scot,P.gamma.f.Scot)
offset.index.Eng<-c(1,no.basis+1,no.basis*2+1,no.basis*3+1,no.basis*4+(no.years-2)*2+1,no.basis*4+(no.years-2)*2+no.basis-3+1)
offset.index.Scot<-c(1,no.basis+1,no.basis*2+1,no.basis*3+1,no.basis*4+(no.years-2)*2+1,no.basis*4+(no.years-2)*2+ncol(A2.reduced.Scot)+1)
P.offset.Eng<-build.penalty(P.offset.list.Eng,index=offset.index.Eng,total.col=ncol(DX.all.Eng))
P.offset.Scot<-build.penalty(P.offset.list.Scot,index=offset.index.Scot,total.col=ncol(DX.all.Scot))

alpha.Eng<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index]),offset=c(offsetm[-omit.index],offsetf[-omit.index]),DXX=DX.all.Eng,P.list=list(list(PD[-(1:8),])),P.index=list(1),exp.index=c(1),lambda=c(-335.4200,341.7743),H=P.offset.Eng,fnval=23000,stepmax=50,gradtol=1e-7,steptol=1e-8)

beta.Eng<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index]),offset=c(offsetm[-omit.index],offsetf[-omit.index]),DXX=DX.all.Eng,P.list=list(list(PD[-(1:8),])),P.index=list(no.basis*2+1),exp.index=c(1),lambda=c(1,5),H=P.offset.Eng,fnval=23000,stepmax=10,gradtol=1e-8,steptol=1e-12)
beta.Eng2<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index]),offset=c(offsetm[-omit.index],offsetf[-omit.index]),DXX=DX.all.Eng,P.list=list(list(PD[-(1:8),])),P.index=list(no.basis*2+1),exp.index=c(1),lambda=c(-7,13),H=P.offset.Eng,fnval=23000,stepmax=5,gradtol=1e-7,steptol=1e-12)

lambda.weightstd<-exp(-50+80*seq(0,1,length=nrow(PD[-(1:8),])))
P.beta.PD.ha<-as.numeric(lambda.weightstd^0.5)*PD[-(1:8),]
P.beta.PD<-build.penalty(list(P.beta.PD.ha),index=no.basis*2+1,total.col=ncol(DX.all.Eng))


all.Eng.no.betaPD<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index]),offset=c(offsetm[-omit.index],offsetf[-omit.index]),DXX=DX.all.Eng,P.list=P.list.Eng[-6],P.index=P.index.Eng[-6],exp.index=c(1,1,1,0,0,0,0),lambda=c(Eng$sp[1:2],Engf$sp[1:2],alpha.Eng$sp,Eng$sp[3],Engf$sp[3],Eng$sp[4],Engf$sp[4]),H=P.beta.PD,fnval=27000,stepmax=2,gradtol=1e-7,steptol=1e-8)
all.Eng.no.betaPD2<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index]),offset=c(offsetm[-omit.index],offsetf[-omit.index]),DXX=DX.all.Eng,P.list=P.list.Eng[-6],P.index=P.index.Eng[-6],exp.index=c(1,1,1,0,0,0,0),lambda=all.Eng.no.betaPD$sp,H=P.beta.PD,fnval=27000,stepmax=2,gradtol=1e-7,steptol=1e-8)
all.Eng.no.betaPD.46<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index]),offset=c(offsetm[-omit.index],offsetf[-omit.index]),DXX=DX.all.Eng,P.list=P.list.Eng[-6],P.index=P.index.Eng[-6],exp.index=c(1,1,1,0,0,0,0),lambda=c(all.Eng.46$sp[-c(9,10)]),H=P.beta.PD,fnval=27000,stepmax=2,gradtol=1e-7,steptol=1e-8)

all.Eng<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index]),offset=c(offsetm[-omit.index],offsetf[-omit.index]),DXX=DX.all.Eng,P.list=P.list.Eng,P.index=P.index.Eng,exp.index=c(1,1,1,0,0,1,0,0),lambda=c(Eng$sp[1:2],Engf$sp[1:2],alpha.Eng$sp,Eng$sp[3],Engf$sp[3],beta.Eng$sp,Eng$sp[4],Engf$sp[4]),H=NULL,typsize=abs(c(Eng$sp[1:2],Engf$sp[1:2],alpha.Eng$sp,Eng$sp[3],Engf$sp[3],beta.Eng$sp,Eng$sp[4],Engf$sp[4])),fscale=27000,stepmax=0.5)
all.Eng2<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index]),offset=c(offsetm[-omit.index],offsetf[-omit.index]),DXX=DX.all.Eng,P.list=P.list.Eng,P.index=P.index.Eng,exp.index=c(1,1,1,0,0,1,0,0),lambda=c(Eng$sp[1:2],Engf$sp[1:2],-7,13,Eng$sp[3],Engf$sp[3],beta.Eng$sp,Eng$sp[4],Engf$sp[4]),H=NULL,typsize=abs(c(Eng$sp[1:2],Engf$sp[1:2],-7,13,Eng$sp[3],Engf$sp[3],beta.Eng$sp,Eng$sp[4],Engf$sp[4])),fscale=40000,stepmax=50)

all.Eng.40<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index]),offset=c(offsetm[-omit.index],offsetf[-omit.index]),DXX=DX.all.Eng,P.list=P.list.Eng,P.index=P.index.Eng,exp.index=c(1,1,1,0,0,1,0,0),lambda=all.Eng$sp,H=NULL,typsize=abs(all.Eng$sp),fscale=40000,stepmax=0.5)
all.Eng.46<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index]),offset=c(offsetm[-omit.index],offsetf[-omit.index]),DXX=DX.all.Eng,P.list=P.list.Eng,P.index=P.index.Eng,exp.index=c(1,1,1,0,0,1,0,0),lambda=all.Eng$sp,H=NULL,typsize=abs(all.Eng$sp),fscale=40000,stepmax=0.5)
all.Eng.46.2<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index]),offset=c(offsetm[-omit.index],offsetf[-omit.index]),DXX=DX.all.Eng,P.list=P.list.Eng,P.index=P.index.Eng,exp.index=c(1,1,1,0,0,1,0,0),lambda=c(all.Eng.46$sp),H=NULL,fnval=23000,stepmax=4,gradtol=1e-8,steptol=1e-12)

all.Eng.56<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index]),offset=c(offsetm[-omit.index],offsetf[-omit.index]),DXX=DX.all.Eng,P.list=P.list.Eng,P.index=P.index.Eng,exp.index=c(1,1,1,0,0,1,0,0),lambda=c(Eng$sp[1:2],Engf$s[1:2],-300,310,Eng$sp[3],Engf$sp[3],-50,80,Eng$sp[4],Engf$sp[4]),H=NULL,fnval=27000,stepmax=2,gradtol=1e-6,steptol=1e-8)
all.Eng.56.3<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index]),offset=c(offsetm[-omit.index],offsetf[-omit.index]),mu=all.Eng.56.2$fv,DXX=DX.all.Eng,P.list=P.list.Eng,P.index=P.index.Eng,exp.index=c(1,1,1,0,0,1,0,0),lambda=c(all.Eng.56.2$sp[1:8],-41,70,all.Eng.56.2$sp[11:12]),H=NULL,fnval=29000,stepmax=10,gradtol=1e-7,steptol=1e-8)

plot(A%*%all.Eng$b[no.basis*2+1:no.basis],type="l")
lines(A%*%all.Eng$b[no.basis*3+1:no.basis],type="l")
lines(A%*%Eng$b[no.basis+1:no.basis],type="l",lty=2)
lines(A%*%Engf$b[no.basis+1:no.basis],type="l",lty=2)

alpha.Scot<-newton.raphson.bic(data=c(datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all.Scot,P.list=list(list(PD[-(1:8),])),P.index=list(1),exp.index=c(1),lambda=c(-10,15),H=P.offset.Scot,fnval=23000,stepmax=50,gradtol=1e-7,steptol=1e-8)

beta.Scot<-newton.raphson.bic(data=c(datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all.Scot,P.list=list(list(PD[-(1:8),])),P.index=list(no.basis*2+1),exp.index=c(1),lambda=c(-2,10),H=P.offset.Scot,typsize=c(5,5),fscale=14000,stepmax=2)
beta.Scot2<-newton.raphson.bic(data=c(datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all.Scot,P.list=list(list(PD[-(1:8),])),P.index=list(no.basis*2+1),exp.index=c(1),lambda=c(-30,50),H=P.offset.Scot,fnval=14513,stepmax=10,gradtol=1e-5,steptol=1e-8)

all.Scot<-newton.raphson.bic(data=c(datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all.Scot,P.list=P.list.Scot,P.index=P.index.Scot,exp.index=c(1,1,1,0,0,1,0,0),lambda=c(Scot$sp[1:2],Scotf$sp[1:2],alpha.Scot2$sp,Scot$sp[3],Scotf$sp[3],beta.Scot$sp,Scot$sp[4],Scotf$sp[4]),H=NULL,typsize=c(0.5,0.5,0.5,0.5,1,1,2,2,1,1,2,2),fscale=1e3,stepmax=50)

lambda.weightstd<-exp(-50+80*seq(0,1,length=nrow(PD[-(1:8),])))
P.beta.PD.ha<-as.numeric(lambda.weightstd^0.5)*PD[-(1:8),]
P.beta.PD.SC<-build.penalty(list(P.beta.PD.ha),index=no.basis*2+1,total.col=ncol(DX.all.Scot))

all.Scot.no.betaPD<-newton.raphson.bic(data=c(datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all.Scot,P.list=P.list.Scot[-6],P.index=P.index.Scot[-6],exp.index=c(1,1,1,0,0,0,0),lambda=c(Scot$sp[1:2],Scotf$sp[1:2],alpha.Scot$sp,Scot$sp[3],Scotf$sp[3],Scot$sp[4],Scotf$sp[4]),H=P.beta.PD.SC,fnval=14000,stepmax=5,gradtol=1e-7,steptol=1e-8)

all.Scot.40<-newton.raphson.bic(data=c(datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all.Scot,P.list=P.list.Scot,P.index=P.index.Scot,exp.index=c(1,1,1,0,0,1,0,0),lambda=c(all.Scot$sp),H=NULL,typsize=abs(all.Scot$sp),fscale=1e3,stepmax=50)
all.Scot.46<-newton.raphson.bic(data=c(datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all.Scot,P.list=P.list.Scot,P.index=P.index.Scot,exp.index=c(1,1,1,0,0,1,0,0),lambda=c(all.Scot$sp),H=NULL,typsize=abs(all.Scot$sp),fscale=1e3,stepmax=50)
all.Scot.56<-newton.raphson.bic(data=c(datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all.Scot,P.list=P.list.Scot,P.index=P.index.Scot,exp.index=c(1,1,1,0,0,1,0,0),lambda=c(Scot$sp[1:2],Scotf$sp[1:2],alpha.Scot$sp,Scot$sp[3],Scotf$sp[3],beta.Scot$sp,Scot$sp[4],Scotf$sp[4]),H=NULL,typsize=abs(c(Scot$sp[1:2],Scotf$sp[1:2],alpha.Scot$sp,Scot$sp[3],Scotf$sp[3],beta.Scot$sp,Scot$sp[4],Scotf$sp[4])),fscale=1e3,stepmax=50)
all.Scot.56.2<-newton.raphson.bic(data=c(datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),mu=all.Scot.56$fv,DXX=DX.all.Scot,P.list=P.list.Scot,P.index=P.index.Scot,exp.index=c(1,1,1,0,0,1,0,0),lambda=all.Scot.56$sp,H=NULL,fnval=10000,stepmax=10,gradtol=1e-7,steptol=1e-8)

######Same sex eng and Scot
m.all.common<-newton.raphson.bic(data=c(datam[-omit.index],datam.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetm.Scot[-omit.index.Scot]),DXX=DX.sex,P.list=P.list.sex,P.index=P.index.sex,exp.index=c(1,0,0),lambda=c(2,4,6,6),typsize=c(1,1,1,1),fscale=1e3,stepmax=30)
f.all.common<-newton.raphson.bic(data=c(dataf[-omit.index],dataf.Scot[-omit.index.Scot]),offset=c(offsetf[-omit.index],offsetf.Scot[-omit.index.Scot]),DXX=DX.sex,P.list=P.list.sex,P.index=P.index.sex,exp.index=c(1,0,0),lambda=c(2,4,6,6),typsize=c(1,1,1,1),fscale=1e3,stepmax=30)

m.all.common.40<-newton.raphson.bic(data=c(datam[-omit.index],datam.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetm.Scot[-omit.index.Scot]),DXX=DX.sex,P.list=P.list.sex,P.index=P.index.sex,exp.index=c(1,0,0),lambda=c(2,4,6,6),typsize=c(1,1,1,1),fscale=1e3,stepmax=30)
f.all.common.40<-newton.raphson.bic(data=c(dataf[-omit.index],dataf.Scot[-omit.index.Scot]),offset=c(offsetf[-omit.index],offsetf.Scot[-omit.index.Scot]),DXX=DX.sex,P.list=P.list.sex,P.index=P.index.sex,exp.index=c(1,0,0),lambda=c(2,4,6,6),typsize=c(1,1,1,1),fscale=1e3,stepmax=30)

lambda.weightsm<-exp(m.all.common$sp[1]+m.all.common$sp[2]*seq(0,1,length=nrow(P)))
lambda.weightsf<-exp(f.all.common$sp[1]+f.all.common$sp[2]*seq(0,1,length=nrow(P)))
P.alpha.m<-as.numeric(lambda.weightsm^0.5)*P
P.alpha.f<-as.numeric(lambda.weightsf^0.5)*P
P.beta.m<-as.numeric(exp(m.all.common$sp[3])^0.5)*P
P.beta.f<-as.numeric(exp(f.all.common$sp[3])^0.5)*P
P.gamma.m<-as.numeric(exp(m.all.common$sp[4])^0.5)*PC
P.gamma.f<-as.numeric(exp(f.all.common$sp[4])^0.5)*PC
P.gamma.m.Scot<-as.numeric(exp(m.all.common$sp[4])^0.5)*PC.Scot
P.gamma.f.Scot<-as.numeric(exp(f.all.common$sp[4])^0.5)*PC.Scot
P.offset.list.m<-list(P.alpha.m,P.alpha.m,P.beta.m,P.beta.m,P.gamma.m,P.gamma.m.Scot)
P.offset.list.f<-list(P.alpha.f,P.alpha.f,P.beta.f,P.beta.f,P.gamma.f,P.gamma.f.Scot)
offset.index<-c(1,no.basis+1,no.basis*2+1,no.basis*3+1,no.basis*4+(no.years-2)*2+1,no.basis*4+(no.years-2)*2+no.basis-3+1)
P.offset.m<-build.penalty(P.offset.list.m,index=offset.index,total.col=ncol(DX.sex))
P.offset.f<-build.penalty(P.offset.list.f,index=offset.index,total.col=ncol(DX.sex))

m.alpha.PD1<-newton.raphson.bic(data=c(datam[-omit.index],datam.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetm.Scot[-omit.index.Scot]),DXX=DX.sex,P.list=list(list(PD1)),P.index=list(1),exp.index=c(0),lambda=c(7),H=P.offset.m,typsize=c(1),fscale=1e3,stepmax=30)
f.alpha.PD1<-newton.raphson.bic(data=c(dataf[-omit.index],dataf.Scot[-omit.index.Scot]),offset=c(offsetf[-omit.index],offsetf.Scot[-omit.index.Scot]),DXX=DX.sex,P.list=list(list(PD1)),P.index=list(1),exp.index=c(0),lambda=c(7),H=P.offset.f,typsize=c(1),fscale=1e3,stepmax=30)
m.beta.PD1<-newton.raphson.bic(data=c(datam[-omit.index],datam.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetm.Scot[-omit.index.Scot]),DXX=DX.sex,P.list=list(list(PD1)),P.index=list(no.basis*2+1),exp.index=c(0),lambda=c(11),H=P.offset.m,typsize=c(1),fscale=1e3,stepmax=2)
f.beta.PD1<-newton.raphson.bic(data=c(dataf[-omit.index],dataf.Scot[-omit.index.Scot]),offset=c(offsetf[-omit.index],offsetf.Scot[-omit.index.Scot]),DXX=DX.sex,P.list=list(list(PD1)),P.index=list(no.basis*2+1),exp.index=c(0),lambda=c(7),H=P.offset.f,typsize=c(1),fscale=1e3,stepmax=2)

all<-newton.raphson.bic(data=c(datam[-omit.index],datam.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetm.Scot[-omit.index.Scot]),DXX=DX.sex,P.list=P.list.PD1,P.index=P.index.PD1,exp.index=c(1,0,0,0,0),lambda=c(m.all.common$sp[1:2],m.alpha.PD1$sp,m.all.common$sp[3],m.beta.PD1$sp,m.all.common$sp[4]),H=NULL,typsize=c(1,1,1,1.5,1,1.5),fscale=1e3,stepmax=3)
all.f<-newton.raphson.bic(data=c(dataf[-omit.index],dataf.Scot[-omit.index.Scot]),offset=c(offsetf[-omit.index],offsetf.Scot[-omit.index.Scot]),DXX=DX.sex,P.list=P.list.PD1,P.index=P.index.PD1,exp.index=c(1,0,0,0,0),lambda=c(f.all.common$sp[1:2],f.alpha.PD1$sp,f.all.common$sp[3],f.beta.PD1$sp,f.all.common$sp[4]),H=NULL,typsize=c(1,1,1,1.5,1,1.5),fscale=1e3,stepmax=3)

all.56<-newton.raphson.bic(data=c(datam[-omit.index],datam.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetm.Scot[-omit.index.Scot]),DXX=DX.sex,P.list=P.list.PD1,P.index=P.index.PD1,exp.index=c(1,0,0,0,0),lambda=c(m.all.common$sp[1:2],m.alpha.PD1$sp,m.all.common$sp[3],m.beta.PD1$sp,m.all.common$sp[4]),H=NULL,typsize=c(1,1,1,1.5,1,1.5),fscale=1e3,stepmax=3)
all.f.56<-newton.raphson.bic(data=c(dataf[-omit.index],dataf.Scot[-omit.index.Scot]),offset=c(offsetf[-omit.index],offsetf.Scot[-omit.index.Scot]),DXX=DX.sex,P.list=P.list.PD1,P.index=P.index.PD1,exp.index=c(1,0,0,0,0),lambda=c(f.all.common$sp[1:2],f.alpha.PD1$sp,f.all.common$sp[3],f.beta.PD1$sp,f.all.common$sp[4]),H=NULL,typsize=c(1,1,1,1.5,1,1.5),fscale=1e3,stepmax=3)

all.40<-newton.raphson.bic(data=c(datam[-omit.index],datam.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetm.Scot[-omit.index.Scot]),DXX=DX.sex,P.list=P.list.PD1,P.index=P.index.PD1,exp.index=c(1,0,0,0,0),lambda=all$sp,H=NULL,typsize=abs(all$sp),fscale=1e3,stepmax=3)
all.f.40<-newton.raphson.bic(data=c(dataf[-omit.index],dataf.Scot[-omit.index.Scot]),offset=c(offsetf[-omit.index],offsetf.Scot[-omit.index.Scot]),DXX=DX.sex,P.list=P.list.PD1,P.index=P.index.PD1,exp.index=c(1,0,0,0,0),lambda=all.f$sp,H=NULL,typsize=abs(all.f$sp),fscale=1e3,stepmax=3)
all.f.40.2<-newton.raphson.bic(data=c(dataf[-omit.index],dataf.Scot[-omit.index.Scot]),offset=c(offsetf[-omit.index],offsetf.Scot[-omit.index.Scot]),DXX=DX.sex,P.list=P.list.PD1,P.index=P.index.PD1,exp.index=c(1,0,0,0,0),lambda=c(all.f$sp[1:4],8,7),H=NULL,typsize=abs(c(all.f$sp[1:4],8,7)),fscale=1e3,stepmax=2)

all.PD.PD1<-newton.raphson.bic(data=c(datam[-omit.index],datam.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetm.Scot[-omit.index.Scot]),DXX=DX.sex,P.list=P.list.PD.PD1,P.index=P.index.PD.PD1,exp.index=c(1,1,0,0,1,0,0),lambda=c(m.all.common$sp[1:2],1,7,m.alpha.PD1$sp,m.all.common$sp[3],1,7,m.beta.PD1$sp,m.all.common$sp[4]),H=NULL,typsize=abs(c(m.all.common$sp[1:2],1,7,m.alpha.PD1$sp,m.all.common$sp[3],1,7,m.beta.PD1$sp,m.all.common$sp[4])),fscale=1e3,stepmax=3)

#######everything tgt
P.index.Eng.Scot<-list(c(1,no.basis*2+1),
			c(no.basis+1,no.basis*3+1),
			c(1,no.basis*2+1),
			c(1,no.basis+1),
			c(no.basis*4+1,no.basis*6+1),
			c(no.basis*5+1,no.basis*7+1),
			c(no.basis*4+1,no.basis*6+1),
			c(no.basis*4+1,no.basis*5+1),
			c(no.basis*8+(no.years-2)*4+1,no.basis*8+(no.years-2)*4+(no.basis-3)*2+1),
			c(no.basis*8+(no.years-2)*4+no.basis-3+1,no.basis*8+(no.years-2)*4+(no.basis-3)*2+ncol(A2.reduced.Scot)+1))
P.list.Eng.Scot=list(list(P,P),list(P,P),list(PD[-(1:8),],PD[-(1:8),]),list(PD1.joint,PD1.joint),list(P,P),list(P,P),list(PD[-(1:8),],PD[-(1:8),]),list(PD1.joint,PD1.joint),list(PC,PC.Scot),list(PC,PC.Scot))

P.index.Eng.Scot2<-list(c(1,no.basis*2+1),
			c(no.basis+1,no.basis*3+1),
			c(1),
			c(no.basis*2+1),
			c(1),
			c(no.basis+1),
			c(no.basis*4+1,no.basis*6+1),
			c(no.basis*5+1,no.basis*7+1),
			c(no.basis*4+1),
			c(no.basis*6+1),
			c(no.basis*4+1),
			c(no.basis*5+1),
			c(no.basis*8+(no.years-2)*4+1,no.basis*8+(no.years-2)*4+(no.basis-3)*2+1),
			c(no.basis*8+(no.years-2)*4+no.basis-3+1,no.basis*8+(no.years-2)*4+(no.basis-3)*2+ncol(A2.reduced.Scot)+1))
P.list.Eng.Scot2=list(list(P,P),list(P,P),list(PD[-(1:8),]),list(PD[-(1:8),]),list(PD1.joint),list(PD1.joint),list(P,P),list(P,P),list(PD[-(1:8),]),list(PD[-(1:8),]),list(PD1.joint),list(PD1.joint),list(PC,PC.Scot),list(PC,PC.Scot))

lambda.weightsm<-exp(m.all.common$sp[1]+m.all.common$sp[2]*seq(0,1,length=nrow(P)))
lambda.weightsf<-exp(f.all.common$sp[1]+f.all.common$sp[2]*seq(0,1,length=nrow(P)))
P.alpha.m<-as.numeric(lambda.weightsm^0.5)*P
P.alpha.f<-as.numeric(lambda.weightsf^0.5)*P
P.beta.m<-as.numeric(exp(m.all.common$sp[3])^0.5)*P
P.beta.f<-as.numeric(exp(f.all.common$sp[3])^0.5)*P
P.gamma.m<-as.numeric(exp(m.all.common$sp[4])^0.5)*PC
P.gamma.f<-as.numeric(exp(f.all.common$sp[4])^0.5)*PC
P.gamma.m.Scot<-as.numeric(exp(m.all.common$sp[4])^0.5)*PC.Scot
P.gamma.f.Scot<-as.numeric(exp(f.all.common$sp[4])^0.5)*PC.Scot
P.offset.list.joint<-list(P.alpha.m,P.alpha.f,P.alpha.m,P.alpha.f,P.beta.m,P.beta.f,P.beta.m,P.beta.f,P.gamma.m,P.gamma.f,P.gamma.m.Scot,P.gamma.f.Scot)
offset.index.joint<-c(1,no.basis+1,no.basis*2+1,no.basis*3+1,no.basis*4+1,no.basis*5+1,no.basis*6+1,no.basis*7+1,no.basis*8+(no.years-2)*4+1,no.basis*8+(no.years-2)*4+no.basis-3+1,no.basis*8+(no.years-2)*4+(no.basis-3)*2+1,no.basis*8+(no.years-2)*4+(no.basis-3)*2+ncol(A2.reduced.Scot)+1)
P.offset.joint<-build.penalty(P.offset.list.joint,index=offset.index.joint,total.col=ncol(DX.all))

j.alpha.PD1<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=list(list(PD1.joint,PD1.joint)),P.index=list(c(1,no.basis+1)),exp.index=c(0),lambda=c(8),H=P.offset.joint,typsize=c(1),fscale=1e3,stepmax=10)
j.alpha.PD1.m<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=list(list(PD1.joint)),P.index=list(c(1)),exp.index=c(0),lambda=c(8),H=P.offset.joint,typsize=c(1),fscale=1e3,stepmax=2)
j.alpha.PD1.f<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=list(list(PD1.joint)),P.index=list(c(no.basis+1)),exp.index=c(0),lambda=c(8),H=P.offset.joint,typsize=c(1),fscale=1e3,stepmax=2)
j.beta.PD1<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=list(list(PD1.joint,PD1.joint)),P.index=list(c(no.basis*4+1,no.basis*5+1)),exp.index=c(0),lambda=c(8),H=P.offset.joint,typsize=c(1),fscale=1e3,stepmax=2)
j.beta.PD1.m<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=list(list(PD1.joint)),P.index=list(c(no.basis*4+1)),exp.index=c(0),lambda=c(8),H=P.offset.joint,typsize=c(1),fscale=1e3,stepmax=2)
j.beta.PD1.f<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=list(list(PD1.joint)),P.index=list(c(no.basis*5+1)),exp.index=c(0),lambda=c(8),H=P.offset.joint,typsize=c(1),fscale=1e3,stepmax=2)
 
j.alpha.PD<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=list(list(PD[-(1:8),],PD[-(1:8),])),P.index=list(c(1,no.basis*2+1)),exp.index=c(1),lambda=c(-7,13),H=P.offset.joint,typsize=c(1,1),fscale=59400,stepmax=50)
j.alpha.PD2<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=list(list(PD[-(1:8),],PD[-(1:8),])),P.index=list(c(1,no.basis*2+1)),exp.index=c(1),lambda=c(2,0),H=P.offset.joint,typsize=c(1,1),fscale=1e3,stepmax=50)
j.alpha.PD.m<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=list(list(PD[-(1:8),])),P.index=list(c(1)),exp.index=c(1),lambda=c(2,0),H=P.offset.joint,typsize=c(1,1),fscale=1e3,stepmax=10)
j.alpha.PD.f<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=list(list(PD[-(1:8),])),P.index=list(c(no.basis*2+1)),exp.index=c(1),lambda=c(-100,120),H=P.offset.joint,typsize=c(1,1),fscale=1e3,stepmax=100)

j.beta.PD<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=list(list(PD[-(1:8),],PD[-(1:8),])),P.index=list(c(no.basis*4+1,no.basis*6+1)),exp.index=c(1),lambda=c(-10,20),H=P.offset.joint,typsize=abs(c(10,30)),fscale=1e3,stepmax=50)
j.beta.PD.ew<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=list(list(PD[-(1:8),])),P.index=list(no.basis*4+1),exp.index=c(1),lambda=c(2,6),H=P.offset.joint,typsize=c(1,1),fscale=1e3,stepmax=5)
j.beta.PD.sc<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=list(list(PD[-(1:8),])),P.index=list(no.basis*6+1),exp.index=c(1),lambda=c(1,8),H=P.offset.joint,typsize=c(1,1),fscale=1e3,stepmax=10)
j.beta.PD.sc.2<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=list(list(PD[-(1:8),])),P.index=list(no.basis*6+1),exp.index=c(1),lambda=c(2,0),H=P.offset.joint,typsize=c(1,1),fscale=1e3,stepmax=10)

P.alpha.PD1<-as.numeric(exp(j.alpha.PD1$sp)^0.5)*PD1.joint
P.beta.PD1<-as.numeric(exp(j.beta.PD1$sp)^0.5)*PD1.joint
P.offset.extra.list.joint<-list(P.alpha.PD1,P.alpha.PD1,P.beta.PD1,P.beta.PD1)
offset.extra.index.joint<-c(1,no.basis+1,no.basis*4+1,no.basis*5+1)
P.offset.extra.joint<-build.penalty(P.offset.extra.list.joint,index=offset.extra.index.joint,total.col=ncol(DX.all))

#lambda.weights.PD<-exp(j.alpha.PD$sp[1]+j.alpha.PD$sp[2]*seq(0,1,length=nrow(PD[-(1:8),])))
#P.alpha.PD<-as.numeric(lambda.weights.PD^0.5)*PD[-(1:8),]
#P.offset.extra.list.joint<-list(P.alpha.PD,P.alpha.PD,P.alpha.PD1,P.alpha.PD1,P.beta.PD1,P.beta.PD1)
#offset.extra.index.joint<-c(1,no.basis*2+1,1,no.basis+1,no.basis*4+1,no.basis*5+1)
#P.offset.extra.joint<-build.penalty(P.offset.extra.list.joint,index=offset.extra.index.joint,total.col=ncol(DX.all))

j.beta.PD<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=list(list(PD[-(1:8),],PD[-(1:8),])),P.index=list(c(no.basis*4+1,no.basis*6+1)),exp.index=c(1),lambda=c(2,5),H=rbind(P.offset.joint,P.offset.extra.joint),typsize=abs(c(0.5,1)),fscale=1e3,stepmax=30)

all.Eng.Scot<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=P.list.Eng.Scot,P.index=P.index.Eng.Scot,exp.index=c(1,1,1,0,0,0,1,0,0,0),lambda=c(m.all.common$sp[1:2],f.all.common$sp[1:2],j.alpha.PD$sp,j.alpha.PD1$sp,m.all.common$sp[3],f.all.common$sp[3],j.beta.PD$sp,j.beta.PD1$sp,m.all.common$sp[4],f.all.common$sp[4]),H=NULL,typsize=abs(c(m.all.common$sp[1:2],f.all.common$sp[1:2],j.alpha.PD$sp,j.alpha.PD1$sp,m.all.common$sp[3],f.all.common$sp[3],j.beta.PD$sp,j.beta.PD1$sp,m.all.common$sp[4],f.all.common$sp[4])),fscale=1,stepmax=10)
all.Eng.Scot<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=P.list.Eng.Scot,P.index=P.index.Eng.Scot,exp.index=c(1,1,1,0,0,0,1,0,0,0),lambda=c(m.all.common$sp[1:2],f.all.common$sp[1:2],-7,13,j.alpha.PD1$sp,m.all.common$sp[3],f.all.common$sp[3],j.beta.PD$sp,j.beta.PD1$sp,m.all.common$sp[4],f.all.common$sp[4]),H=NULL,typsize=abs(c(m.all.common$sp[1:2],f.all.common$sp[1:2],7,13,j.alpha.PD1$sp,m.all.common$sp[3],f.all.common$sp[3],j.beta.PD$sp,j.beta.PD1$sp,m.all.common$sp[4],f.all.common$sp[4])),fscale=59460,stepmax=10)
all.Eng.Scot.50<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=P.list.Eng.Scot,P.index=P.index.Eng.Scot,exp.index=c(1,1,1,0,0,0,1,0,0,0),lambda=all.Eng.Scot$sp,H=NULL,typsize=abs(all.Eng.Scot$sp),fscale=44000,stepmax=10)

all.Eng.Scot.40<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=P.list.Eng.Scot,P.index=P.index.Eng.Scot,exp.index=c(1,1,1,0,0,0,1,0,0,0),lambda=c(all.Eng.Scot$sp[1:9],-7,13,all.Eng.Scot$sp[12:14]),H=NULL,typsize=abs(all.Eng.Scot$sp),fscale=59460,stepmax=10)
all.Eng.Scot.52<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=P.list.Eng.Scot,P.index=P.index.Eng.Scot,exp.index=c(1,1,1,0,0,0,1,0,0,0),lambda=all.Eng.Scot$sp,H=NULL,typsize=abs(all.Eng.Scot$sp),fscale=44000,stepmax=10)

all.Eng.Scot.56<-newton.raphson.bic(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]),DXX=DX.all,P.list=P.list.Eng.Scot,P.index=P.index.Eng.Scot,exp.index=c(1,1,1,0,0,0,1,0,0,0),lambda=c(m.all.common$sp[1:2],f.all.common$sp[1:2],j.alpha.PD$sp,j.alpha.PD1$sp,m.all.common$sp[3],f.all.common$sp[3],j.beta.PD$sp,j.beta.PD1$sp,m.all.common$sp[4],f.all.common$sp[4]),H=NULL,typsize=abs(c(m.all.common$sp[1:2],f.all.common$sp[1:2],j.alpha.PD$sp,j.alpha.PD1$sp,m.all.common$sp[3],f.all.common$sp[3],j.beta.PD$sp,j.beta.PD1$sp,m.all.common$sp[4],f.all.common$sp[4])),fscale=44000,stepmax=10)

save(m.all.common,f.all.common,m.alpha.PD1,m.beta.PD1,f.alpha.PD1,f.beta.PD1,all,all.f,fit,fit.f,fit.m.sep,fit.m.Scot.sep,fit.f.sep,fit.f.Scot.sep,file="same sp just PD1 pivot.RData")
save(Eng,Scot,Engf,Scotf,alpha.Eng,alpha.Scot,beta.Eng,beta.Scot,all.Eng,all.Scot,fit.Eng,fit.Scot,fit.m.sep,fit.m.Scot.sep,fit.f.sep,fit.f.Scot.sep,file="males and females Eng and Scot pivot.RData")
save(m.all.common,f.all.common,j.alpha.PD,j.alpha.PD1,j.beta.PD,j.beta.PD1,all.Eng.Scot,fit.Eng.Scot,fit.m.sep,fit.m.Scot.sep,fit.f.sep,fit.f.Scot.sep,file="same sp Eng Scot all tgt pivot.RData")
save(Eng.40,Scot.40,Engf.40,Scotf.40,all.Eng.40,all.Scot.40,file="males and females Eng and Scot pivot 40 years.RData")
save(Eng,Scot,Engf,Scotf,all.56,all.f.56,fit,fit.f,fit.m.sep,fit.m.Scot.sep,fit.f.sep,fit.f.Scot.sep,file="males and females Eng and Scot pivot 56 years.RData")

#####males and females
par(mar=c(5,5,4,2),mfrow=c(1,2))
plot(A%*%all.Eng.no.betaPD$b[no.basis+1:no.basis],type="l",col=2,cex.lab=1.5,cex.main=1.5,main="Basline Mortality Schedule (EW)",xlab="Age",ylab=bquote(s[alpha](x)))
lines(A%*%all.Eng.no.betaPD$b[1:no.basis],type="l",col=4)
lines(A%*%Eng$b[1:no.basis],type="l",col=4,lty=2)
lines(A%*%Engf$b[1:no.basis],type="l",col=2,lty=2)
legend("topleft",c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,2,4,2),x.intersp=0.3,y.intersp=0.7)

plot(A%*%all.Eng.no.betaPD$b[no.basis*2+1:no.basis],type="l",col=4,cex.lab=1.5,cex.main=1.5,main="Age-specific Improvement Rates (EW)",xlab="Age",ylab=bquote(s[beta](x)))
lines(A%*%all.Eng.no.betaPD$b[no.basis*3+1:no.basis],type="l",col=2)
lines(A%*%Eng$b[no.basis+1:no.basis],type="l",col=4,lty=2)
lines(A%*%Engf$b[no.basis+1:no.basis],type="l",col=2,lty=2)
legend(30,-1.3,c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,2,4,2),x.intersp=0.3,y.intersp=0.7)

plot(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%all.Eng.no.betaPD$b[no.basis*4+1:(no.years-2)],type="l",col=4,cex.lab=1.5,cex.main=1.5,main="Period Effect (EW)",xlab="Year",ylab=bquote(kappa[t]))
lines(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%all.Eng.no.betaPD$b[no.basis*4+no.years-2+1:(no.years-2)],type="l",col=2)
lines(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%Eng$b[no.basis*2+1:(no.years-2)],type="l",col=4,lty=2)
lines(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%Engf$b[no.basis*2+1:(no.years-2)],type="l",col=2,lty=2)
legend(1965,-0.04,c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,2,4,2),x.intersp=0.3,y.intersp=0.7)

plot((1961-last.age):(1961+no.years-2),c(rep(0,cohort.omit),A2%*%(C.Q.trans%*%all.Eng.no.betaPD$b[no.basis*4+(no.years-2)*2+1:(no.basis-3)])),type="l",col=4,cex.lab=1.5,cex.main=1.5,main="Cohort Effect (EW)",xlab="Year of Birth",ylab=bquote(s[gamma](t-x)))
lines((1961-last.age):(1961+no.years-2),c(rep(0,cohort.omit),A2%*%(C.Q.trans%*%all.Eng.no.betaPD$b[no.basis*4+(no.years-2)*2+no.basis-3+1:(no.basis-3)])),type="l",col=2)
lines((1961-last.age):(1961+no.years-2),c(rep(0,cohort.omit),A2%*%(C.Q.trans%*%Eng$b[no.basis*2+(no.years-2)+1:(no.basis-3)])),type="l",col=4,lty=2)
lines((1961-last.age):(1961+no.years-2),c(rep(0,cohort.omit),A2%*%(C.Q.trans%*%Engf$b[no.basis*2+(no.years-2)+1:(no.basis-3)])),type="l",col=2,lty=2)
legend(1840,-0.1,c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,2,4,2),x.intersp=0.3,y.intersp=0.7)


par(mar=c(5,5,4,2),mfrow=c(1,2))
plot(A.Scot%*%all.Scot$b[no.basis+1:no.basis],type="l",col=3,cex.lab=1.5,cex.main=1.5,main="Basline Mortality Schedule (Scotland)",xlab="Age",ylab=bquote(s[alpha](x)))
lines(A.Scot%*%all.Scot$b[1:no.basis],type="l",col=1)
lines(A.Scot%*%Scot$b[1:no.basis],type="l",col=1,lty=2)
lines(A.Scot%*%Scotf$b[1:no.basis],type="l",col=3,lty=2)
legend("topleft",c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(1,3,1,3),x.intersp=0.3,y.intersp=0.7)

plot(A.Scot%*%all.Scot$b[no.basis*2+1:no.basis],type="l",col=1,cex.lab=1.5,cex.main=1.5,main="Age-specific Improvement Rates (Scotland)",xlab="Age",ylab=bquote(s[beta](x)),ylim=c(-2.3,0.2))
lines(A.Scot%*%all.Scot$b[no.basis*3+1:no.basis],type="l",col=3)
lines(A.Scot%*%Scot$b[no.basis+1:no.basis],type="l",col=1,lty=2)
lines(A.Scot%*%Scotf$b[no.basis+1:no.basis],type="l",col=3,lty=2)
legend(30,-1.3,c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(1,3,1,3),x.intersp=0.3,y.intersp=0.7)

plot(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%all.Scot$b[no.basis*4+1:(no.years-2)],type="l",col=1,cex.lab=1.5,cex.main=1.5,main="Period Effect (Scotland)",xlab="Year",ylab=bquote(kappa[t]))
lines(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%all.Scot$b[no.basis*4+no.years-2+1:(no.years-2)],type="l",col=3)
lines(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%Scot$b[no.basis*2+1:(no.years-2)],type="l",col=1,lty=2)
lines(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%Scotf$b[no.basis*2+1:(no.years-2)],type="l",col=3,lty=2)
legend(1965,-0.05,c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(1,3,1,3),x.intersp=0.3,y.intersp=0.7)

plot((1961-last.age.Scot):(1961+no.years-2),c(rep(0,cohort.omit),A2.Scot%*%(C.Q.trans.Scot%*%all.Scot$b[no.basis*4+(no.years-2)*2+1:(ncol(A2.Scot)-3)])),type="l",col=1,cex.lab=1.5,cex.main=1.5,main="Cohort Effect (Scotland)",xlab="Year of Birth",ylab=bquote(s[gamma](t-x)))
lines((1961-last.age.Scot):(1961+no.years-2),c(rep(0,cohort.omit),A2.Scot%*%(C.Q.trans.Scot%*%all.Scot$b[no.basis*4+(no.years-2)*2+ncol(A2.Scot)-3+1:(ncol(A2.Scot)-3)])),type="l",col=3)
lines((1961-last.age.Scot):(1961+no.years-2),c(rep(0,cohort.omit),A2.Scot%*%(C.Q.trans.Scot%*%Scot$b[no.basis*2+(no.years-2)+1:(ncol(A2.Scot)-3)])),type="l",col=1,lty=2)
lines((1961-last.age.Scot):(1961+no.years-2),c(rep(0,cohort.omit),A2.Scot%*%(C.Q.trans.Scot%*%Scotf$b[no.basis*2+(no.years-2)+1:(ncol(A2.Scot)-3)])),type="l",col=3,lty=2)
legend(1840,-0.1,c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(1,3,1,3),x.intersp=0.3,y.intersp=0.7)

#####Eng and Scot
par(mar=c(5,5,4,2),mfrow=c(1,2))
plot(A%*%all.56$b[no.basis+1:no.basis],type="l",col=1,cex.lab=1.5,cex.main=1.5,main="Basline Mortality Schedule (Males)",xlab="Age",ylab=bquote(s[alpha](x)))
lines(A%*%all.56$b[1:no.basis],type="l",col=4)
lines(A%*%Eng$b[1:no.basis],type="l",col=4,lty=2)
lines(A%*%Scot$b[1:no.basis],type="l",col=1,lty=2)
legend("topleft",c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,1,4,1),x.intersp=0.3,y.intersp=0.7)

plot(A%*%all.56$b[no.basis*3+1:no.basis],type="l",ylim=c(-2.2,0.3),col=1,cex.lab=1.5,cex.main=1.5,main="Age-specific Improvement Rates (Males)",xlab="Age",ylab=bquote(s[beta](x)))
lines(A%*%all.56$b[no.basis*2+1:no.basis],type="l",col=4)
lines(A%*%Eng$b[no.basis+1:no.basis],type="l",col=4,lty=2)
lines(A%*%Scot$b[no.basis+1:no.basis],type="l",col=1,lty=2)
legend(30,-1.3,c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,1,4,1),x.intersp=0.3,y.intersp=0.7)

plot(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%all.56$b[no.basis*4+1:(no.years-2)],type="l",col=4,cex.lab=1.5,cex.main=1.5,main="Period Effect (Males)",xlab="Year",ylab=bquote(kappa[t]))
lines(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%all.56$b[no.basis*4+no.years-2+1:(no.years-2)],type="l",col=1)
lines(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%Eng$b[no.basis*2+1:(no.years-2)],type="l",col=4,lty=2)
lines(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%Scot$b[no.basis*2+1:(no.years-2)],type="l",col=1,lty=2)
legend(1965,-0.04,c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,1,4,1),x.intersp=0.3,y.intersp=0.7)

plot((1961-last.age):(1961+no.years-2),c(rep(NA,last.age-last.age.Scot),rep(0,cohort.omit),A2.Scot%*%(C.Q.trans.Scot%*%all.56$b[no.basis*4+(no.years-2)*2+no.basis-3+1:(ncol(A2.Scot)-3)])),type="l",col=1,ylim=c(-0.35,0.25),cex.lab=1.5,cex.main=1.5,main="Cohort Effect (Males)",xlab="Year of Birth",ylab=bquote(s[gamma](t-x)))
lines((1961-last.age):(1961+no.years-2),c(rep(0,cohort.omit),A2%*%(C.Q.trans%*%all.56$b[no.basis*4+(no.years-2)*2+1:(no.basis-3)])),type="l",col=4)
lines((1961-last.age):(1961+no.years-2),c(rep(0,cohort.omit),A2%*%(C.Q.trans%*%Eng$b[no.basis*2+(no.years-2)+1:(no.basis-3)])),type="l",col=4,lty=2)
lines((1961-last.age):(1961+no.years-2),c(rep(NA,last.age-last.age.Scot),rep(0,cohort.omit),A2.Scot%*%(C.Q.trans.Scot%*%Scot$b[no.basis*2+(no.years-2)+1:(ncol(A2.Scot)-3)])),type="l",col=1,lty=2)
legend(1840,-0.1,c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,1,4,1),x.intersp=0.3,y.intersp=0.7)


par(mar=c(5,5,4,2),mfrow=c(1,2))
plot(A%*%all.f.56$b[no.basis+1:no.basis],type="l",col=6,cex.lab=1.5,cex.main=1.5,main="Basline Mortality Schedule (Females)",xlab="Age",ylab=bquote(s[alpha](x)))
lines(A%*%all.f.56$b[1:no.basis],type="l",col=2)
lines(A%*%Engf$b[1:no.basis],type="l",col=2,lty=2)
lines(A%*%Scotf$b[1:no.basis],type="l",col=6,lty=2)
legend("topleft",c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(2,6,2,6),x.intersp=0.3,y.intersp=0.7)

plot(A%*%all.f.56$b[no.basis*3+1:no.basis],type="l",ylim=c(-2.2,0.3),col=6,cex.lab=1.5,cex.main=1.5,main="Age-specific Improvement Rates (Females)",xlab="Age",ylab=bquote(s[beta](x)))
lines(A%*%all.f.56$b[no.basis*2+1:no.basis],type="l",col=2)
lines(A%*%Engf$b[no.basis+1:no.basis],type="l",col=2,lty=2)
lines(A%*%Scotf$b[no.basis+1:no.basis],type="l",col=6,lty=2)
legend(30,-1.3,c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(2,6,2,6),x.intersp=0.3,y.intersp=0.7)

plot(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%all.f.56$b[no.basis*4+1:(no.years-2)],type="l",col=2,cex.lab=1.5,cex.main=1.5,main="Period Effect (Females)",xlab="Year",ylab=bquote(kappa[t]))
lines(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%all.f.56$b[no.basis*4+no.years-2+1:(no.years-2)],type="l",col=6)
lines(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%Engf$b[no.basis*2+1:(no.years-2)],type="l",col=2,lty=2)
lines(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%Scotf$b[no.basis*2+1:(no.years-2)],type="l",col=6,lty=2)
legend(1965,-0.04,c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(2,6,2,6),x.intersp=0.3,y.intersp=0.7)

plot((1961-last.age):(1961+no.years-2),c(rep(0,cohort.omit),A2%*%(C.Q.trans%*%all.f.56$b[no.basis*4+(no.years-2)*2+1:(no.basis-3)])),type="l",ylim=c(-0.13,0.16),col=2,cex.lab=1.5,cex.main=1.5,main="Cohort Effect (Females)",xlab="Year of Birth",ylab=bquote(s[gamma](t-x)))
lines((1961-last.age):(1961+no.years-2),c(rep(NA,last.age-last.age.Scot),rep(0,cohort.omit),A2.Scot%*%(C.Q.trans.Scot%*%all.f.56$b[no.basis*4+(no.years-2)*2+no.basis-3+1:(ncol(A2.Scot)-3)])),type="l",col=6)
lines((1961-last.age):(1961+no.years-2),c(rep(0,cohort.omit),A2%*%(C.Q.trans%*%Engf$b[no.basis*2+(no.years-2)+1:(no.basis-3)])),type="l",col=2,lty=2)
lines((1961-last.age):(1961+no.years-2),c(rep(NA,last.age-last.age.Scot),rep(0,cohort.omit),A2.Scot%*%(C.Q.trans.Scot%*%Scotf$b[no.basis*2+(no.years-2)+1:(ncol(A2.Scot)-3)])),type="l",col=6,lty=2)
legend(1840,-0.1,c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(2,6,2,6),x.intersp=0.3,y.intersp=0.7)


#####all tgt
par(mar=c(5,5,4,2),mfrow=c(1,4))
plot(A%*%all.Eng.Scot$b[1:no.basis],type="l",col=4,ylim=c(-9,0),cex.lab=1.5,cex.main=1.5,main="Basline Mortality Schedule \n(EW Males)",xlab="Age",ylab=bquote(s[alpha](x)))
lines(A%*%Eng$b[1:no.basis],type="l",col=4,lty=2)

plot(A%*%all.Eng.Scot$b[no.basis*2+1:no.basis],type="l",col=1,cex.lab=1.5,cex.main=1.5,main="Basline Mortality Schedule \n(Scotland Males)",xlab="Age",ylab=bquote(s[alpha](x)))
lines(A%*%Scot$b[1:no.basis],type="l",col=1,lty=2)

plot(A%*%all.Eng.Scot$b[no.basis+1:no.basis],type="l",col=2,cex.lab=1.5,cex.main=1.5,main="Basline Mortality Schedule \n(EW Females)",xlab="Age",ylab=bquote(s[alpha](x)))
lines(A%*%Engf$b[1:no.basis],type="l",col=2,lty=2)

plot(A%*%all.Eng.Scot$b[no.basis*3+1:no.basis],type="l",col=6,cex.lab=1.5,cex.main=1.5,main="Basline Mortality Schedule \n(Scotland Females)",xlab="Age",ylab=bquote(s[alpha](x)))
lines(A%*%Scotf$b[1:no.basis],type="l",col=6,lty=2)
legend(23,-8,c("Joint","Separate"),bty="n",lty=c(1,2),col=c(1,1),x.intersp=0.3,seg.len=1)

par(mar=c(5,5,4,2),mfrow=c(1,4))
plot(A%*%all.Eng.Scot$b[no.basis*4+1:no.basis],type="l",col=4,cex.lab=1.5,cex.main=1.5,main="Age-specific Improvement Rates \n(EW Males)",xlab="Age",ylab=bquote(s[beta](x)))
lines(A%*%Eng$b[no.basis+1:no.basis],type="l",col=4,lty=2)

plot(A%*%Scot$b[no.basis+1:no.basis],type="l",lty=2,col=1,cex.lab=1.5,cex.main=1.5,main="Age-specific Improvement Rates \n(Scotland Males)",xlab="Age",ylab=bquote(s[beta](x)))
lines(A%*%all.Eng.Scot$b[no.basis*6+1:no.basis],type="l",col=1,lty=1)

plot(A%*%all.Eng.Scot$b[no.basis*5+1:no.basis],type="l",col=2,cex.lab=1.5,cex.main=1.5,main="Age-specific Improvement Rates \n(EW Females)",xlab="Age",ylab=bquote(s[beta](x)))
lines(A%*%Engf$b[no.basis+1:no.basis],type="l",col=2,lty=2)

plot(A%*%all.Eng.Scot$b[no.basis*7+1:no.basis],type="l",col=6,cex.lab=1.5,cex.main=1.5,main="Age-specific Improvement Rates \n(Scotland Females)",xlab="Age",ylab=bquote(s[beta](x)))
lines(A%*%Scotf$b[no.basis+1:no.basis],type="l",col=6,lty=2)
legend(23,-1.5,c("Joint","Separate"),bty="n",lty=c(1,2),col=c(1,1),x.intersp=0.3,seg.len=1)

par(mar=c(5,5,4,2),mfrow=c(1,4))
plot(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%all.Eng.Scot$b[no.basis*8+1:(no.years-2)],type="l",col=4,cex.lab=1.5,cex.main=1.5,main="Period Effect (EW Males)",xlab="Year",ylab=bquote(kappa[t]))
lines(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%Eng$b[no.basis*2+1:(no.years-2)],type="l",col=4,lty=2)

plot(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%all.Eng.Scot$b[no.basis*8+(no.years-2)*2+1:(no.years-2)],type="l",col=1,cex.lab=1.5,cex.main=1.5,main="Period Effect (Scotland Males)",xlab="Year",ylab=bquote(kappa[t]))
lines(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%Scot$b[no.basis*2+1:(no.years-2)],type="l",col=1,lty=2)

plot(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%all.Eng.Scot$b[no.basis*8+no.years-2+1:(no.years-2)],type="l",col=2,cex.lab=1.5,cex.main=1.5,main="Period Effect (EW Females)",xlab="Year",ylab=bquote(kappa[t]))
lines(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%Engf$b[no.basis*2+1:(no.years-2)],type="l",col=2,lty=2)

plot(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%all.Eng.Scot$b[no.basis*8+(no.years-2)*3+1:(no.years-2)],type="l",col=6,cex.lab=1.5,cex.main=1.5,main="Period Effect (Scotland Females)",xlab="Year",ylab=bquote(kappa[t]))
lines(1961:(1961+no.years-1),inv.kappa.trans[,-c(1,2)]%*%Scotf$b[no.basis*2+1:(no.years-2)],type="l",col=6,lty=2)
legend(1970,-0.055,c("Joint","Separate"),bty="n",lty=c(1,2),col=c(1,1),x.intersp=0.3,seg.len=1)

par(mar=c(5,5,4,2),mfrow=c(1,4))
plot((1961-last.age):(1961+no.years-2),c(rep(0,cohort.omit),A2%*%(C.Q.trans%*%all.Eng.Scot$b[no.basis*8+(no.years-2)*4+1:(no.basis-3)])),type="l",col=4,cex.lab=1.5,cex.main=1.5,main="Cohort Effect (EW Males)",xlab="Year of Birth",ylab=bquote(s[gamma](t-x)))
lines((1961-last.age):(1961+no.years-2),c(rep(0,cohort.omit),A2%*%(C.Q.trans%*%Eng$b[no.basis*2+(no.years-2)+1:(no.basis-3)])),type="l",col=4,lty=2)

plot((1961-last.age):(1961+no.years-2),c(rep(NA,last.age-last.age.Scot),rep(0,cohort.omit),A2.Scot%*%(C.Q.trans.Scot%*%Scot$b[no.basis*2+(no.years-2)+1:(ncol(A2.Scot)-3)])),type="l",lty=2,col=1,cex.lab=1.5,cex.main=1.5,main="Cohort Effect (Scotland Males)",xlab="Year of Birth",ylab=bquote(s[gamma](t-x)))
lines((1961-last.age):(1961+no.years-2),c(rep(NA,last.age-last.age.Scot),rep(0,cohort.omit),A2.Scot%*%(C.Q.trans.Scot%*%all.Eng.Scot$b[no.basis*8+(no.years-2)*4+(no.basis-3)*2+1:(ncol(A2.Scot)-3)])),type="l",col=1,lty=1)

plot((1961-last.age):(1961+no.years-2),c(rep(0,cohort.omit),A2%*%(C.Q.trans%*%Engf$b[no.basis*2+(no.years-2)+1:(no.basis-3)])),type="l",lty=2,col=2,cex.lab=1.5,cex.main=1.5,main="Cohort Effect (EW Females)",xlab="Year of Birth",ylab=bquote(s[gamma](t-x)))
lines((1961-last.age):(1961+no.years-2),c(rep(0,cohort.omit),A2%*%(C.Q.trans%*%all.Eng.Scot$b[no.basis*8+(no.years-2)*4+no.basis-3+1:(no.basis-3)])),type="l",col=2,lty=1)

plot((1961-last.age):(1961+no.years-2),c(rep(NA,last.age-last.age.Scot),rep(0,cohort.omit),A2.Scot%*%(C.Q.trans.Scot%*%all.Eng.Scot$b[no.basis*8+(no.years-2)*4+(no.basis-3)*2+ncol(A2.Scot)-3+1:(ncol(A2.Scot)-3)])),type="l",lty=1,col=6,cex.lab=1.5,cex.main=1.5,main="Cohort Effect (Scotland Females)",xlab="Year of Birth",ylab=bquote(s[gamma](t-x)))
lines((1961-last.age):(1961+no.years-2),c(rep(NA,last.age-last.age.Scot),rep(0,cohort.omit),A2.Scot%*%(C.Q.trans.Scot%*%Scotf$b[no.basis*2+(no.years-2)+1:(ncol(A2.Scot)-3)])),type="l",col=6,lty=2)
legend(1830,-0.13,c("Joint","Separate"),bty="n",lty=c(1,2),col=c(1,1),x.intersp=0.3,seg.len=1)


IRLS<-function (data,mu=data,DXX,B,offset) {
mu<-ifelse(mu==0,1e-5,mu)
dev<-ll.sat<-sum(dpois(data,data,log=TRUE))
eta<-log(mu)

converged<-FALSE
while(!converged) {
	z<-(data-mu)/mu+eta-offset
	w<-as.numeric(mu^0.5)
	#fit<-lm(c(w*z,rep(0,nrow(B)))~rbind(w*DXX,B)-1)	

X.QR<-qr(w*DXX)
QR.pivot<-X.QR$pivot
Q<-qr.Q(X.QR)
R<-qr.R(X.QR)
total.col<-ncol(DXX)
SVD<-svd(rbind(R,B[,QR.pivot]))
#####################lowest svd d
svd.index<-which(SVD$d > max(SVD$d)*sqrt(.Machine$double.eps))
#####################
#svd.index<-which(SVD$d > max(SVD$d)*1e-7)
V<-SVD$v[,svd.index]
U1<-SVD$u[1:ncol(DXX),svd.index]
D<-SVD$d[svd.index]
coef<-(V%*%diag(1/D)%*%t(U1)%*%(t(Q)%*%(w*z)))[sort.list(QR.pivot)]
	eta<-DXX%*%coef+offset
#	eta<-DXX%*%coef(fit)+offset
	mu<-exp(eta)
	old.dev<-dev
	dev<-2*(ll.sat-sum(dpois(data,mu,log=TRUE)))
	if(abs(dev-old.dev)<1e-6*old.dev) converged<-TRUE
}
b=coef

#X.QR<-qr(w*DXX)
#QR.pivot<-X.QR$pivot
#Q<-qr.Q(X.QR)
#R<-qr.R(X.QR)
#total.col<-ncol(DXX)
#SVD<-svd(rbind(R,B[,QR.pivot]))
#####################lowest svd d
#svd.index<-which(SVD$d > max(SVD$d)*sqrt(.Machine$double.eps))
#####################
#svd.index<-which(SVD$d > max(SVD$d)*1e-7)
#U1<-SVD$u[1:ncol(DXX),svd.index]

#b=coef(fit)
fv=mu
edf<-sum(diag(crossprod(U1)))
rss<-sum((w*(z-DXX%*%coef))^2)
#rss<-sum(w*(z-DXX%*%coef(fit))^2)
bic<-dev+log(length(data))*edf
bic2<-rss+log(length(data))*edf
cred.int<-crossprod((1/D)*t(V))[sort.list(QR.pivot),sort.list(QR.pivot)]
list(b=b,fv=fv,dev=dev,bic=bic,bic2=bic2,cred.int=cred.int)
}

build.penalty<-function(P,index,total.col){
B<-c()
for(i in 1:length(P)) {
B<-rbind(B,cbind(matrix(0,nrow(P[[i]]),index[i]-1),P[[i]],matrix(0,nrow(P[[i]]),total.col-ncol(P[[i]])-index[i]+1)))
	}
B
}

last.age=104;no.years=56;start.year=1961;end.year=start.year+no.years-1
last.age.Scot=99;no.years.Scot=no.years;start.year.Scot=1961;end.year.Scot=start.year.Scot+no.years.Scot-1
cohort.omit=6

no.basis=40;extra.knots=0
knots<-seq(0,1,length=no.basis-2)
dk<-knots[2]-knots[1]	
knots<-c(knots[1]-dk*(3:1),knots,knots[no.basis-2]+dk*(1:(3+extra.knots)))
age<-seq(0,1,length=last.age)
age.Scot<-((1:last.age.Scot)-1)/(last.age-1)
cohort.years<-seq(0,1,length=last.age+no.years-1-cohort.omit)
cohort.years.Scot<-seq(0,1,length=last.age.Scot+no.years.Scot-1-cohort.omit)

bspline<-function (x,k,i,m=2) {
if (m==-1) {basis<-as.numeric(x<k[i+1] & x>=k[i])} else {
	z0<-(x-k[i])/(k[i+m+1]-k[i])
	z1<-(k[i+m+2]-x)/(k[i+m+2]-k[i+1])
	basis<-z0*bspline(x,k,i,m-1)+z1*bspline(x,k,i+1,m-1) }
basis
}

timeindex<-((1:no.years)-1)/(no.years-1)
timeindex<-timeindex-mean(timeindex)
timeindex.Scot<-timeindex[(start.year.Scot-start.year.Scot+1):(end.year.Scot-start.year.Scot+1)]

A<-c()
for(j in 1:(no.basis+extra.knots)) {
A<-cbind(A,bspline(age,knots,j))
}

A.Scot<-c()
for(j in 1:(no.basis+extra.knots)) {
A.Scot<-cbind(A.Scot,bspline(age.Scot,knots,j))
}

#Scot.index<-which(apply(A.Scot,2,sum)==0)
#A.Scot<-A.Scot[,-Scot.index]

DX<-rep(1,no.years)%x%A
DX.Scot<-rep(1,no.years.Scot)%x%A.Scot

DX.period<-timeindex%x%A
DX.period.Scot<-timeindex.Scot%x%A.Scot

kappa.trans<-rbind(rep(1,no.years),timeindex,cbind(matrix(0,no.years-2,2),diag(no.years-2)))
inv.kappa.trans<-solve(kappa.trans)
kappa.mat<-rbind(t(solve(kappa.trans)[1,3:no.years])%x%rep(1,last.age),t(solve(kappa.trans)[2,3:no.years])%x%rep(1,last.age),diag(no.years-2)%x%rep(1,last.age))

kappa.trans.Scot<-rbind(rep(1,no.years.Scot),timeindex.Scot,cbind(matrix(0,no.years.Scot-2,2),diag(no.years.Scot-2)))
inv.kappa.trans.Scot<-solve(kappa.trans.Scot)
kappa.mat.Scot<-rbind(t(solve(kappa.trans.Scot)[1,3:no.years.Scot])%x%rep(1,last.age.Scot),t(solve(kappa.trans.Scot)[2,3:no.years.Scot])%x%rep(1,last.age.Scot),diag(no.years.Scot-2)%x%rep(1,last.age.Scot))

A2<-c()
for(j in 1:no.basis) {
A2<-cbind(A2,bspline(cohort.years,knots,j))
}

CX.constraints<-rbind(A2[1,],apply(A2,2,sum),A2[nrow(A2),])
C.Q<-qr.Q(qr(t(CX.constraints)),complete=T)
C.Q.trans<-C.Q[,-(1:3)]
A2.reduced<-A2%*%C.Q.trans
A2.reduced<-rbind(matrix(0,cohort.omit,no.basis-3),A2.reduced)

A2.Scot<-A2[((start.year.Scot-last.age.Scot)-(start.year-last.age)+1):((start.year.Scot-last.age.Scot)-(start.year-last.age)+no.years.Scot+last.age.Scot-1-cohort.omit),]
A2.Scot<-A2.Scot[,colSums(A2.Scot)!=0]
CX.constraints.Scot<-rbind(A2.Scot[1,],apply(A2.Scot,2,sum),A2.Scot[nrow(A2.Scot),])
C.Q.Scot<-qr.Q(qr(t(CX.constraints.Scot)),complete=T)
C.Q.trans.Scot<-C.Q.Scot[,-(1:3)]
A2.reduced.Scot<-A2.Scot%*%C.Q.trans.Scot
A2.reduced.Scot<-rbind(matrix(0,cohort.omit,ncol(A2.Scot)-3),A2.reduced.Scot)

##########
A2.SC<-c()
for(j in 1:no.basis) {
A2.SC<-cbind(A2.SC,bspline(cohort.years.Scot,knots,j))
}

CX.constraints.SC<-rbind(A2.SC[1,],apply(A2.SC,2,sum),A2.SC[nrow(A2.SC),])
C.Q.SC<-qr.Q(qr(t(CX.constraints.SC)),complete=T)
C.Q.trans.SC<-C.Q.SC[,-(1:3)]
A2.reduced.SC<-A2.SC%*%C.Q.trans.SC
A2.reduced.SC<-rbind(matrix(0,cohort.omit,no.basis-3),A2.reduced.SC)
###########

CX<-c()
for(i in 1:no.years){
CX<-rbind(CX,A2.reduced[(i-1)+(last.age:1),])
}

CX2<-c()
for(i in 1:no.years.Scot){
CX2<-rbind(CX2,A2.reduced.Scot[(i-1)+(last.age.Scot:1),])
}

DX.cohort<-CX
DX.cohort.Scot<-CX2

CX2.SC<-c()
for(i in 1:no.years.Scot){
CX2.SC<-rbind(CX2.SC,A2.reduced.SC[(i-1)+(last.age.Scot:1),])
}

DX.cohort.SC<-CX2.SC

omit.index<-c();for(i in 1:cohort.omit){omit.index<-c(omit.index,last.age*(i-1)+(last.age-cohort.omit+i):last.age)}
omit.index.Scot<-c();for(i in 1:cohort.omit){omit.index.Scot<-c(omit.index.Scot,last.age.Scot*(i-1)+(last.age.Scot-cohort.omit+i):last.age.Scot)}

DX.all<-cbind(adiag(DX[-omit.index,],DX[-omit.index,],DX.Scot[-omit.index.Scot,],DX.Scot[-omit.index.Scot,]),adiag(DX.period[-omit.index,],DX.period[-omit.index,],DX.period.Scot[-omit.index.Scot,],DX.period.Scot[-omit.index.Scot,]),adiag(kappa.mat[-omit.index,],kappa.mat[-omit.index,],kappa.mat.Scot[-omit.index.Scot,],kappa.mat.Scot[-omit.index.Scot,]),adiag(DX.cohort[-omit.index,],DX.cohort[-omit.index,],DX.cohort.Scot[-omit.index.Scot,],DX.cohort.Scot[-omit.index.Scot,]))
DX.sex<-cbind(adiag(DX[-omit.index,],DX.Scot[-omit.index.Scot,]),adiag(DX.period[-omit.index,],DX.period.Scot[-omit.index.Scot,]),adiag(kappa.mat[-omit.index,],kappa.mat.Scot[-omit.index.Scot,]),adiag(DX.cohort[-omit.index,],DX.cohort.Scot[-omit.index.Scot,]))
DX.all.Eng<-cbind(adiag(DX[-omit.index,],DX[-omit.index,]),adiag(DX.period[-omit.index,],DX.period[-omit.index,]),adiag(kappa.mat[-omit.index,],kappa.mat[-omit.index,]),adiag(DX.cohort[-omit.index,],DX.cohort[-omit.index,]))
DX.all.Scot<-cbind(adiag(DX.Scot[-omit.index.Scot,],DX.Scot[-omit.index.Scot,]),adiag(DX.period.Scot[-omit.index.Scot,],DX.period.Scot[-omit.index.Scot,]),adiag(kappa.mat.Scot[-omit.index.Scot,],kappa.mat.Scot[-omit.index.Scot,]),adiag(DX.cohort.Scot[-omit.index.Scot,],DX.cohort.Scot[-omit.index.Scot,]))
####cohort with no.basis
DX.all.SC<-cbind(adiag(DX.Scot[-omit.index.Scot,],DX.Scot[-omit.index.Scot,]),adiag(DX.period.Scot[-omit.index.Scot,],DX.period.Scot[-omit.index.Scot,]),adiag(kappa.mat.Scot[-omit.index.Scot,],kappa.mat.Scot[-omit.index.Scot,]),adiag(DX.cohort.SC[-omit.index.Scot,],DX.cohort.SC[-omit.index.Scot,]))
####

DX.all.single<-cbind(DX,DX.period,kappa.mat,DX.cohort)[-omit.index,]
DX.all.single.Scot<-cbind(DX.Scot,DX.period.Scot,kappa.mat.Scot,DX.cohort.Scot)[-omit.index.Scot,]

P<-diff(diag(no.basis+extra.knots),differences=2)
PD<-t(c(1,-1))%x%diag(no.basis+extra.knots)
PD1<-t(c(1,-1))%x%diff(diag(no.basis+extra.knots),differences=1)
PD1.joint<-t(c(1,0,-1))%x%diff(diag(no.basis+extra.knots),differences=1)
PC<-diff(diag(no.basis),differences=2)%*%C.Q.trans
PC.Scot<-diff(diag(ncol(A2.Scot)),differences=2)%*%C.Q.trans.Scot
PC.SC<-diff(diag(ncol(A2.SC)),differences=2)%*%C.Q.trans.SC

P.index.PD1<-list(c(1,no.basis+extra.knots+1),
			c(1),
			c((no.basis+extra.knots)*2+1,(no.basis+extra.knots)*3+1),
			c((no.basis+extra.knots)*2+1),
			c((no.basis+extra.knots)*4+(no.years-2)*2+1,(no.basis+extra.knots)*4+(no.years-2)*2+no.basis-3+1))

P.index.Eng<-list(1,
			no.basis+extra.knots+1,
			1,
			(no.basis+extra.knots)*2+1,
			(no.basis+extra.knots)*3+1,
			(no.basis+extra.knots)*2+1,
			(no.basis+extra.knots)*4+(no.years-2)*2+1,
			(no.basis+extra.knots)*4+(no.years-2)*2+no.basis-3+1)
P.index.Scot<-list(1,
			no.basis+extra.knots+1,
			1,
			(no.basis+extra.knots)*2+1,
			(no.basis+extra.knots)*3+1,
			(no.basis+extra.knots)*2+1,
			(no.basis+extra.knots)*4+(no.years-2)*2+1,
			(no.basis+extra.knots)*4+(no.years-2)*2+ncol(A2.reduced.Scot)+1)
P.index.SC<-list(1,
			no.basis+extra.knots+1,
			1,
			(no.basis+extra.knots)*2+1,
			(no.basis+extra.knots)*3+1,
			(no.basis+extra.knots)*2+1,
			(no.basis+extra.knots)*4+(no.years-2)*2+1,
			(no.basis+extra.knots)*4+(no.years-2)*2+no.basis-3+1)
P.index.Eng.Scot<-list(c(1,(no.basis+extra.knots)*2+1),
			c((no.basis+extra.knots)+1,(no.basis+extra.knots)*3+1),
			c(1,(no.basis+extra.knots)*2+1),
			c(1,(no.basis+extra.knots)+1),
			c((no.basis+extra.knots)*4+1,(no.basis+extra.knots)*6+1),
			c((no.basis+extra.knots)*5+1,(no.basis+extra.knots)*7+1),
			c((no.basis+extra.knots)*4+1,(no.basis+extra.knots)*6+1),
			c((no.basis+extra.knots)*4+1,(no.basis+extra.knots)*5+1),
			c((no.basis+extra.knots)*8+(no.years-2)*4+1,(no.basis+extra.knots)*8+(no.years-2)*4+(no.basis-3)*2+1),
			c((no.basis+extra.knots)*8+(no.years-2)*4+no.basis-3+1,(no.basis+extra.knots)*8+(no.years-2)*4+(no.basis-3)*2+ncol(A2.reduced.Scot)+1))


############Eng and Scot male and female
#lambda<-all.Eng.56.2$sp
lambda<-all.Eng.no.betaPD$sp
lambda<-c(lambda[1:8],-50,80,lambda[9:10])
weights.m<-c(exp(lambda[1]+lambda[2]*seq(0,1,length=no.basis-2)),rep(exp(lambda[1]+lambda[2]),extra.knots))
weights.f<-c(exp(lambda[3]+lambda[4]*seq(0,1,length=no.basis-2)),rep(exp(lambda[3]+lambda[4]),extra.knots))
weights.d<-c(exp(lambda[5]+lambda[6]*seq(0,1,length=no.basis-8)),rep(exp(lambda[5]+lambda[6]),extra.knots))
weights.td<-c(exp(lambda[9]+lambda[10]*seq(0,1,length=no.basis-8)),rep(exp(lambda[9]+lambda[10]),extra.knots))
P.alpha.m<-as.numeric(weights.m^0.5)*P
P.alpha.f<-as.numeric(weights.f^0.5)*P
P.alpha.d<-as.numeric(weights.d^0.5)*PD[-(1:8),]
P.beta.m<-as.numeric(exp(lambda[7])^0.5)*P
P.beta.f<-as.numeric(exp(lambda[8])^0.5)*P
P.beta.d<-as.numeric(weights.td^0.5)*PD[-(1:8),]
P.gamma.m<-as.numeric(exp(lambda[11])^0.5)*PC
P.gamma.f<-as.numeric(exp(lambda[12])^0.5)*PC
P.weighted<-list(P.alpha.m,P.alpha.f,P.alpha.d,P.beta.m,P.beta.f,P.beta.d,P.gamma.m,P.gamma.f)
B<-build.penalty(P.weighted,index=unlist(P.index.Eng),total.col=ncol(DX.all.Eng))

fit.Eng<-IRLS(data=c(datam[-omit.index],dataf[-omit.index]),DXX=DX.all.Eng,B=B,offset=c(offsetm[-omit.index],offsetf[-omit.index]))

lambda<-all.Scot.56$sp
weights.m<-c(exp(lambda[1]+lambda[2]*seq(0,1,length=no.basis-2)),rep(exp(lambda[1]+lambda[2]),extra.knots))
weights.f<-c(exp(lambda[3]+lambda[4]*seq(0,1,length=no.basis-2)),rep(exp(lambda[3]+lambda[4]),extra.knots))
weights.d<-c(exp(lambda[5]+lambda[6]*seq(0,1,length=no.basis-8)),rep(exp(lambda[5]+lambda[6]),extra.knots))
weights.td<-c(exp(lambda[9]+lambda[10]*seq(0,1,length=no.basis-8)),rep(exp(lambda[9]+lambda[10]),extra.knots))
P.alpha.m<-as.numeric(weights.m^0.5)*P
P.alpha.f<-as.numeric(weights.f^0.5)*P
P.alpha.d<-as.numeric(weights.d^0.5)*PD[-(1:8),]
P.beta.m<-as.numeric(exp(lambda[7])^0.5)*P
P.beta.f<-as.numeric(exp(lambda[8])^0.5)*P
P.beta.d<-as.numeric(weights.td^0.5)*PD[-(1:8),]
P.gamma.m<-as.numeric(exp(lambda[11])^0.5)*PC.Scot
P.gamma.f<-as.numeric(exp(lambda[12])^0.5)*PC.Scot
P.gamma.m.SC<-as.numeric(exp(lambda[11])^0.5)*PC.SC
P.gamma.f.SC<-as.numeric(exp(lambda[12])^0.5)*PC.SC
P.weighted<-list(P.alpha.m,P.alpha.f,P.alpha.d,P.beta.m,P.beta.f,P.beta.d,P.gamma.m,P.gamma.f)
B<-build.penalty(P.weighted,index=unlist(P.index.Scot),total.col=ncol(DX.all.Scot))
P.weighted.SC<-list(P.alpha.m,P.alpha.f,P.alpha.d,P.beta.m,P.beta.f,P.beta.d,P.gamma.m.SC,P.gamma.f.SC)
B.SC<-build.penalty(P.weighted.SC,index=unlist(P.index.SC),total.col=ncol(DX.all.SC))

fit.Scot<-IRLS(data=c(datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),DXX=DX.all.Scot,B=B,offset=c(offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]))
fit.SC<-IRLS(data=c(datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),DXX=DX.all.SC,B=B.SC,offset=c(offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]))

##################Eng and Scot same sex
lambda<-all.56$sp
weights.m<-c(exp(lambda[1]+lambda[2]*seq(0,1,length=no.basis-2)),rep(exp(lambda[1]+lambda[2]),extra.knots))
P.alpha.m<-as.numeric(weights.m^0.5)*P
P.alpha.d<-as.numeric(exp(lambda[3])^0.5)*PD1
P.beta.m<-as.numeric(exp(lambda[4])^0.5)*P
P.beta.d<-as.numeric(exp(lambda[5])^0.5)*PD1
P.gamma.m<-as.numeric(exp(lambda[6])^0.5)*PC
P.gamma.m.Scot<-as.numeric(exp(lambda[6])^0.5)*PC.Scot
P.weighted<-list(P.alpha.m,P.alpha.m,P.alpha.d,P.beta.m,P.beta.m,P.beta.d,P.gamma.m,P.gamma.m.Scot)
B<-build.penalty(P.weighted,index=unlist(P.index.PD1),total.col=ncol(DX.sex))

fit<-IRLS(data=c(datam[-omit.index],datam.Scot[-omit.index.Scot]),DXX=DX.sex,B=B,offset=c(offsetm[-omit.index],offsetm.Scot[-omit.index.Scot]))

lambda<-all.f.56$sp
weights.m<-c(exp(lambda[1]+lambda[2]*seq(0,1,length=no.basis-2)),rep(exp(lambda[1]+lambda[2]),extra.knots))
P.alpha.m<-as.numeric(weights.m^0.5)*P
P.alpha.d<-as.numeric(exp(lambda[3])^0.5)*PD1
P.beta.m<-as.numeric(exp(lambda[4])^0.5)*P
P.beta.d<-as.numeric(exp(lambda[5])^0.5)*PD1
P.gamma.m<-as.numeric(exp(lambda[6])^0.5)*PC
P.gamma.m.Scot<-as.numeric(exp(lambda[6])^0.5)*PC.Scot
P.weighted<-list(P.alpha.m,P.alpha.m,P.alpha.d,P.beta.m,P.beta.m,P.beta.d,P.gamma.m,P.gamma.m.Scot)
B<-build.penalty(P.weighted,index=unlist(P.index.PD1),total.col=ncol(DX.sex))

fit.f<-IRLS(data=c(dataf[-omit.index],dataf.Scot[-omit.index.Scot]),DXX=DX.sex,B=B,offset=c(offsetf[-omit.index],offsetf.Scot[-omit.index.Scot]))

lambda.m<-Eng$sp
lambda.f<-Scot$sp
weights.m<-c(exp(lambda.m[1]+lambda.m[2]*seq(0,1,length=no.basis-2)),rep(exp(lambda.m[1]+lambda.m[2]),extra.knots))
weights.f<-c(exp(lambda.f[1]+lambda.f[2]*seq(0,1,length=no.basis-2)),rep(exp(lambda.f[1]+lambda.f[2]),extra.knots))
P.alpha.m<-as.numeric(weights.m^0.5)*P
P.alpha.f<-as.numeric(weights.f^0.5)*P
P.beta.m<-as.numeric(exp(lambda.m[3])^0.5)*P
P.beta.f<-as.numeric(exp(lambda.f[3])^0.5)*P
P.gamma.m<-as.numeric(exp(lambda.m[4])^0.5)*PC
P.gamma.f<-as.numeric(exp(lambda.f[4])^0.5)*PC.Scot
P.weighted.m<-list(P.alpha.m=P.alpha.m,P.beta.m=P.beta.m,P.gamma.m=P.gamma.m)
P.weighted.f<-list(P.alpha.f=P.alpha.f,P.beta.f=P.beta.f,P.gamma.f=P.gamma.f)
B.m<-build.penalty(P.weighted.m,index=c(1,no.basis+extra.knots+1,(no.basis+extra.knots)*2+(no.years-2)+1),total.col=ncol(DX.all.single))
B.f<-build.penalty(P.weighted.f,index=c(1,no.basis+extra.knots+1,(no.basis+extra.knots)*2+(no.years.Scot-2)+1),total.col=ncol(DX.all.single.Scot))

fit.m.sep<-IRLS(data=c(datam[-omit.index]),DXX=DX.all.single,B=B.m,offset=c(offsetm[-omit.index]))
fit.m.Scot.sep<-IRLS(data=c(datam.Scot[-omit.index.Scot]),DXX=DX.all.single.Scot,B=B.f,offset=c(offsetm.Scot[-omit.index.Scot]))

lambda.m<-Engf$sp
lambda.f<-Scotf$sp
weights.m<-c(exp(lambda.m[1]+lambda.m[2]*seq(0,1,length=no.basis-2)),rep(exp(lambda.m[1]+lambda.m[2]),extra.knots))
weights.f<-c(exp(lambda.f[1]+lambda.f[2]*seq(0,1,length=no.basis-2)),rep(exp(lambda.f[1]+lambda.f[2]),extra.knots))
P.alpha.m<-as.numeric(weights.m^0.5)*P
P.alpha.f<-as.numeric(weights.f^0.5)*P
P.beta.m<-as.numeric(exp(lambda.m[3])^0.5)*P
P.beta.f<-as.numeric(exp(lambda.f[3])^0.5)*P
P.gamma.m<-as.numeric(exp(lambda.m[4])^0.5)*PC
P.gamma.f<-as.numeric(exp(lambda.f[4])^0.5)*PC.Scot
P.weighted.m<-list(P.alpha.m=P.alpha.m,P.beta.m=P.beta.m,P.gamma.m=P.gamma.m)
P.weighted.f<-list(P.alpha.f=P.alpha.f,P.beta.f=P.beta.f,P.gamma.f=P.gamma.f)
B.m<-build.penalty(P.weighted.m,index=c(1,no.basis+extra.knots+1,(no.basis+extra.knots)*2+(no.years-2)+1),total.col=ncol(DX.all.single))
B.f<-build.penalty(P.weighted.f,index=c(1,no.basis+extra.knots+1,(no.basis+extra.knots)*2+(no.years.Scot-2)+1),total.col=ncol(DX.all.single.Scot))

fit.f.sep<-IRLS(data=c(dataf[-omit.index]),DXX=DX.all.single,B=B.m,offset=c(offsetf[-omit.index]))
fit.f.Scot.sep<-IRLS(data=c(dataf.Scot[-omit.index.Scot]),DXX=DX.all.single.Scot,B=B.f,offset=c(offsetf.Scot[-omit.index.Scot]))

####all tgt
lambda<-all.Eng.Scot.56.3.alt$sp
lambda<-c(all.Eng.Scot.56.no.betaPD$sp[1:9],-56.03874,86,all.Eng.Scot.56.no.betaPD$sp[10:12])
weights.m<-c(exp(lambda[1]+lambda[2]*seq(0,1,length=no.basis-2)),rep(exp(lambda[1]+lambda[2]),extra.knots))
weights.f<-c(exp(lambda[3]+lambda[4]*seq(0,1,length=no.basis-2)),rep(exp(lambda[3]+lambda[4]),extra.knots))
weights.d<-c(exp(lambda[5]+lambda[6]*seq(0,1,length=no.basis-8)),rep(exp(lambda[5]+lambda[6]),extra.knots))
weights.td<-c(exp(lambda[10]+lambda[11]*seq(0,1,length=no.basis-8)),rep(exp(lambda[10]+lambda[11]),extra.knots))
P.alpha.m<-as.numeric(weights.m^0.5)*P
P.alpha.f<-as.numeric(weights.f^0.5)*P
P.alpha.d<-as.numeric(weights.d^0.5)*PD[-(1:8),]
P.alpha.d1<-as.numeric(exp(lambda[7])^0.5)*PD1.joint
P.beta.m<-as.numeric(exp(lambda[8])^0.5)*P
P.beta.f<-as.numeric(exp(lambda[9])^0.5)*P
P.beta.d<-as.numeric(weights.td^0.5)*PD[-(1:8),]
P.beta.d1<-as.numeric(exp(lambda[12])^0.5)*PD1.joint
P.gamma.m<-as.numeric(exp(lambda[13])^0.5)*PC
P.gamma.f<-as.numeric(exp(lambda[14])^0.5)*PC
P.gamma.m.Scot<-as.numeric(exp(lambda[13])^0.5)*PC.Scot
P.gamma.f.Scot<-as.numeric(exp(lambda[14])^0.5)*PC.Scot
P.weighted<-list(P.alpha.m,P.alpha.m,P.alpha.f,P.alpha.f,P.alpha.d,P.alpha.d,P.alpha.d1,P.alpha.d1,P.beta.m,P.beta.m,P.beta.f,P.beta.f,P.beta.d,P.beta.d,P.beta.d1,P.beta.d1,P.gamma.m,P.gamma.m.Scot,P.gamma.f,P.gamma.f.Scot)
B<-build.penalty(P.weighted,index=unlist(P.index.Eng.Scot),total.col=ncol(DX.all))

fit.Eng.Scot<-IRLS(data=c(datam[-omit.index],dataf[-omit.index],datam.Scot[-omit.index.Scot],dataf.Scot[-omit.index.Scot]),DXX=DX.all,B=B,offset=c(offsetm[-omit.index],offsetf[-omit.index],offsetm.Scot[-omit.index.Scot],offsetf.Scot[-omit.index.Scot]))

A.extra<-c()
for(j in 1:(no.basis+extra.knots)) {
A.extra<-cbind(A.extra,bspline(0:119/(last.age-1),knots,j))
}

#####Males and females joint
par(mar=c(5,5,4,2))
plot(A.extra[,]%*%fit.Eng$b[no.basis+extra.knots+1:(no.basis+extra.knots)],type="l",col=2,cex.lab=1.5,cex.main=1.5,main="Basline Mortality Schedule Extrapolated to 120 (EW)",xlab="Age",ylab=bquote(s[alpha](x)))
lines(A.extra[,]%*%fit.Eng$b[1:(no.basis+extra.knots)],type="l",col=4)
lines(A.extra[,]%*%fit.Eng$b[no.basis+extra.knots+1:(no.basis+extra.knots)],type="l",col=2)
lines(A.extra[,]%*%fit.m.sep$b[1:(no.basis+extra.knots)],type="l",col=4,lty=2)
lines(A.extra[,]%*%fit.f.sep$b[1:(no.basis+extra.knots)],type="l",col=2,lty=2)

lines(A.extra[,]%*%fit.Eng$b[1:(no.basis+extra.knots)]+1.96*sqrt(diag(A.extra%*%fit.Eng$cred.int[1:(no.basis+extra.knots),1:(no.basis+extra.knots)]%*%t(A.extra))),lty=2,col=4)
lines(A.extra[,]%*%fit.Eng$b[1:(no.basis+extra.knots)]-1.96*sqrt(diag(A.extra%*%fit.Eng$cred.int[1:(no.basis+extra.knots),1:(no.basis+extra.knots)]%*%t(A.extra))),lty=2,col=4)
lines(A.extra[,]%*%fit.Eng$b[no.basis+extra.knots+1:(no.basis+extra.knots)]+1.96*sqrt(diag(A.extra%*%fit.Eng$cred.int[no.basis+extra.knots+1:(no.basis+extra.knots),no.basis+extra.knots+1:(no.basis+extra.knots)]%*%t(A.extra))),lty=2,col=2)
lines(A.extra[,]%*%fit.Eng$b[no.basis+extra.knots+1:(no.basis+extra.knots)]-1.96*sqrt(diag(A.extra%*%fit.Eng$cred.int[no.basis+extra.knots+1:(no.basis+extra.knots),no.basis+extra.knots+1:(no.basis+extra.knots)]%*%t(A.extra))),lty=2,col=2)

legend("topleft",c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,2,4,2))

plot(A.extra%*%fit.Eng$b[(no.basis+extra.knots)*2+1:(no.basis+extra.knots)],type="l",col=4,ylim=c(-2,0.5),cex.lab=1.5,cex.main=1.5,main="Age-specific Improvement Rates Extrapolated to 120 (EW)",xlab="Age",ylab=bquote(s[beta](x)))
lines(A.extra%*%fit.Eng$b[(no.basis+extra.knots)*3+1:(no.basis+extra.knots)],type="l",col=2)
lines(A.extra%*%fit.m.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)],type="l",lty=2,col=4)
lines(A.extra%*%fit.f.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)],type="l",col=2,lty=2)

lines(A.extra[,]%*%fit.Eng$b[(no.basis+extra.knots)*2+1:(no.basis+extra.knots)]+1.96*sqrt(diag(A.extra%*%fit.Eng$cred.int[(no.basis+extra.knots)*2+1:(no.basis+extra.knots),(no.basis+extra.knots)*2+1:(no.basis+extra.knots)]%*%t(A.extra))),lty=3,col=4)
lines(A.extra[,]%*%fit.Eng$b[(no.basis+extra.knots)*2+1:(no.basis+extra.knots)]-1.96*sqrt(diag(A.extra%*%fit.Eng$cred.int[(no.basis+extra.knots)*2+1:(no.basis+extra.knots),(no.basis+extra.knots)*2+1:(no.basis+extra.knots)]%*%t(A.extra))),lty=3,col=4)
lines(A.extra[,]%*%fit.Eng$b[(no.basis+extra.knots)*3+1:(no.basis+extra.knots)]+1.96*sqrt(diag(A.extra%*%fit.Eng$cred.int[(no.basis+extra.knots)*3+1:(no.basis+extra.knots),(no.basis+extra.knots)*3+1:(no.basis+extra.knots)]%*%t(A.extra))),lty=3,col=2)
lines(A.extra[,]%*%fit.Eng$b[(no.basis+extra.knots)*3+1:(no.basis+extra.knots)]-1.96*sqrt(diag(A.extra%*%fit.Eng$cred.int[(no.basis+extra.knots)*3+1:(no.basis+extra.knots),(no.basis+extra.knots)*3+1:(no.basis+extra.knots)]%*%t(A.extra))),lty=3,col=2)

lines(A.extra[,]%*%fit.m.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)]+1.96*sqrt(diag(A.extra%*%fit.m.sep$cred.int[(no.basis+extra.knots)+1:(no.basis+extra.knots),(no.basis+extra.knots)+1:(no.basis+extra.knots)]%*%t(A.extra))),lty=4,col=4)
lines(A.extra[,]%*%fit.m.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)]-1.96*sqrt(diag(A.extra%*%fit.m.sep$cred.int[(no.basis+extra.knots)+1:(no.basis+extra.knots),(no.basis+extra.knots)+1:(no.basis+extra.knots)]%*%t(A.extra))),lty=4,col=4)
lines(A.extra[,]%*%fit.f.sep$b[(no.basis+extr]a.knots)+1:(no.basis+extra.knots)]+1.96*sqrt(diag(A.extra%*%fit.f.sep$cred.int[(no.basis+extra.knots)+1:(no.basis+extra.knots),(no.basis+extra.knots)+1:(no.basis+extra.knots)]%*%t(A.extra))),lty=4,col=2)
lines(A.extra[,]%*%fit.f.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)]-1.96*sqrt(diag(A.extra%*%fit.f.sep$cred.int[(no.basis+extra.knots)+1:(no.basis+extra.knots),(no.basis+extra.knots)+1:(no.basis+extra.knots)]%*%t(A.extra))),lty=4,col=2)

legend("topleft",c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,2,4,2))

plot(A.extra[,]%*%fit.Scot$b[no.basis+extra.knots+1:(no.basis+extra.knots)],type="l",col=3,cex.lab=1.5,cex.main=1.5,main="Basline Mortality Schedule Extrapolated to 120 (Scotland)",xlab="Age",ylab=bquote(s[alpha](x)))
lines(A.extra[,]%*%fit.Scot$b[1:(no.basis+extra.knots)],type="l",col=1)
lines(A.extra[,]%*%fit.Scot$b[no.basis+extra.knots+1:(no.basis+extra.knots)],type="l",col=3)
lines(A.extra[,]%*%fit.m.Scot.sep$b[1:(no.basis+extra.knots)],type="l",col=1,lty=2)
lines(A.extra[,]%*%fit.f.Scot.sep$b[1:(no.basis+extra.knots)],type="l",col=3,lty=2)
legend("topleft",c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(1,3,1,3))

plot(A.extra%*%fit.Scot$b[(no.basis+extra.knots)*2+1:(no.basis+extra.knots)],type="l",col=1,ylim=c(-2.3,0.5),cex.lab=1.5,cex.main=1.5,main="Age-specific Improvement Rates Extrapolated to 120 (Scotland)",xlab="Age",ylab=bquote(s[beta](x)))
lines(A.extra%*%fit.Scot$b[(no.basis+extra.knots)*3+1:(no.basis+extra.knots)],type="l",col=3)
lines(A.extra%*%fit.m.Scot.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)],type="l",lty=2,col=1)
lines(A.extra%*%fit.f.Scot.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)],type="l",col=3,lty=2)
legend(85,-1,c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(1,3,1,3))

par(mfrow=c(1,2),mar=c(5,6,4,2),cex.lab=1.8,cex.axis=1.3,cex.main=1.5)
plot(1962:2016,diff(inv.kappa.trans[,-c(1,2)]%*%fit.Eng$b[(no.basis+extra.knots)*4+1:(no.years-2)]),type="l",col=4,xlab="Year",ylab=as.expression(bquote(kappa[t]^"*")),main="Period Effect on difference scale (England and Wales)")
lines(1962:2016,diff(inv.kappa.trans[,-c(1,2)]%*%fit.Eng$b[(no.basis+extra.knots)*4+no.years-2+1:(no.years-2)]),type="l",col=2)
lines(1962:2016,diff(inv.kappa.trans[,-c(1,2)]%*%fit.m.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)]),type="l",col=4,lty=2)
lines(1962:2016,diff(inv.kappa.trans[,-c(1,2)]%*%fit.f.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)]),type="l",col=2,lty=2)
legend(1965,0.067,c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,2,4,2),x.intersp=0.3,y.intersp=0.6)

plot(1864:2015,diff(A2.reduced%*%fit.Eng$b[(no.basis+extra.knots)*4+(no.years-2)*2+1:(no.basis-3)])[-(1:6)],type="l",col=4,xlab="Year of Birth",ylab=as.expression(bquote(s[gamma]^"*"~"(t-x)")),main="Cohort Effect on difference scale (England and Wales)")
lines(1864:2015,diff(A2.reduced%*%fit.Eng$b[(no.basis+extra.knots)*4+(no.years-2)*2+no.basis-3+1:ncol(A2.reduced)])[-(1:6)],type="l",col=2)
lines(1864:2015,diff(A2.reduced%*%fit.m.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(no.basis-3)])[-(1:6)],type="l",col=4,lty=2)
lines(1864:2015,diff(A2.reduced%*%fit.f.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:ncol(A2.reduced)])[-(1:6)],type="l",col=2,lty=2)
legend(1878,0.0170,c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,2,4,2),x.intersp=0.3,y.intersp=0.6)


plot(1857:2015,A2.reduced%*%fit.Eng$b[(no.basis+extra.knots)*4+(no.years-2)*2+1:(no.basis-3)],type="l",col=4,xlab="Year of Birth",ylab=as.expression(bquote(s[gamma]^"*"~"(t-x)")),main="Cohort Effect on difference scale (England and Wales)")
lines(1857:2015,A2.reduced%*%fit.Eng$b[(no.basis+extra.knots)*4+(no.years-2)*2+no.basis-3+1:ncol(A2.reduced)],type="l",col=2)
lines(1857:2015,A2.reduced%*%fit.m.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(no.basis-3)],type="l",col=4,lty=2)
lines(1857:2015,A2.reduced%*%fit.f.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:ncol(A2.reduced)],type="l",col=2,lty=2)
legend(1878,0.0175,c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,2,4,2),x.intersp=0.3,y.intersp=0.6)

par(mfrow=c(1,2),mar=c(5,6,4,2),cex.lab=1.8,cex.axis=1.3,cex.main=1.5)
plot(1962:2016,diff(inv.kappa.trans[,-c(1,2)]%*%fit.Scot$b[(no.basis+extra.knots)*4+1:(no.years-2)]),type="l",col=1,xlab="Year",ylab=as.expression(bquote(kappa[t]^"*")),main="Period Effect on difference scale (Scotland)")
lines(1962:2016,diff(inv.kappa.trans[,-c(1,2)]%*%fit.Scot$b[(no.basis+extra.knots)*4+no.years-2+1:(no.years-2)]),type="l",col=3)
lines(1962:2016,diff(inv.kappa.trans[,-c(1,2)]%*%fit.m.Scot.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)]),type="l",col=1,lty=2)
lines(1962:2016,diff(inv.kappa.trans[,-c(1,2)]%*%fit.f.Scot.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)]),type="l",col=3,lty=2)
legend(1972,0.067,c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(1,3,1,3),x.intersp=0.3,y.intersp=0.6)

plot(1869:2015,diff(A2.reduced.Scot%*%fit.Scot$b[(no.basis+extra.knots)*4+(no.years-2)*2+1:ncol(A2.reduced.Scot)])[-(1:6)],type="l",col=1,xlab="Year of Birth",ylab=as.expression(bquote(s[gamma]^"*"~"(t-x)")),main="Cohort Effect on difference scale (Scotland)")
lines(1869:2015,diff(A2.reduced.Scot%*%fit.Scot$b[(no.basis+extra.knots)*4+(no.years-2)*2+ncol(A2.reduced.Scot)+1:ncol(A2.reduced.Scot)])[-(1:6)],type="l",col=3)
lines(1869:2015,diff(A2.reduced.Scot%*%fit.m.Scot.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:ncol(A2.reduced.Scot)])[-(1:6)],type="l",col=1,lty=2)
lines(1869:2015,diff(A2.reduced.Scot%*%fit.f.Scot.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:ncol(A2.reduced.Scot)])[-(1:6)],type="l",col=3,lty=2)
legend(1850,0.04,c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(1,3,1,3),x.intersp=0.3,y.intersp=0.6)

###England and Scotland joint
par(mar=c(5,5,4,2))
plot(A.extra[,]%*%fit$b[no.basis+extra.knots+1:(no.basis+extra.knots)],type="l",col=1,cex.lab=1.5,cex.main=1.5,main="Basline Mortality Schedule Extrapolated to 120 (Males)",xlab="Age",ylab=bquote(s[alpha](x)))
lines(A.extra[,]%*%fit$b[1:(no.basis+extra.knots)],type="l",col=4)
lines(A.extra[,]%*%fit.m.sep$b[1:(no.basis+extra.knots)],type="l",col=4,lty=2)
lines(A.extra[,]%*%fit.m.Scot.sep$b[1:(no.basis+extra.knots)],type="l",col=1,lty=2)
legend("topleft",c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,1,4,1))

plot(A.extra%*%fit$b[(no.basis+extra.knots)*2+1:(no.basis+extra.knots)],type="l",col=4,ylim=c(-2,0.5),cex.lab=1.5,cex.main=1.5,main="Age-specific Improvement Rates Extrapolated to 120 (Males)",xlab="Age",ylab=bquote(s[beta](x)))
lines(A.extra%*%fit$b[(no.basis+extra.knots)*3+1:(no.basis+extra.knots)],type="l",col=1)
lines(A.extra%*%fit.m.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)],type="l",lty=2,col=4)
lines(A.extra%*%fit.m.Scot.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)],type="l",col=1,lty=2)
legend("topleft",c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,1,4,1))

plot(A.extra[,]%*%fit.f$b[no.basis+extra.knots+1:(no.basis+extra.knots)],type="l",col=3,cex.lab=1.5,cex.main=1.5,main="Basline Mortality Schedule Extrapolated to 120 (Females)",xlab="Age",ylab=bquote(s[alpha](x)))
lines(A.extra[,]%*%fit.f$b[1:(no.basis+extra.knots)],type="l",col=2)
lines(A.extra[,]%*%fit.f.sep$b[1:(no.basis+extra.knots)],type="l",col=2,lty=2)
lines(A.extra[,]%*%fit.f.Scot.sep$b[1:(no.basis+extra.knots)],type="l",col=3,lty=2)
legend("topleft",c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(2,3,2,3))

plot(A.extra%*%fit.f$b[(no.basis+extra.knots)*2+1:(no.basis+extra.knots)],type="l",col=2,ylim=c(-2,0.5),cex.lab=1.5,cex.main=1.5,main="Age-specific Improvement Rates Extrapolated to 120 (Scotland)",xlab="Age",ylab=bquote(s[beta](x)))
lines(A.extra%*%fit.f$b[(no.basis+extra.knots)*3+1:(no.basis+extra.knots)],type="l",col=3)
lines(A.extra%*%fit.f.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)],type="l",lty=2,col=2)
lines(A.extra%*%fit.f.Scot.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)],type="l",col=3,lty=2)
legend("topleft",c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(2,3,2,3))

par(mfrow=c(1,2),cex.axis=1.3,cex.lab=1.5,cex.main=1.5)
plot(1961:2016,inv.kappa.trans[,-c(1,2)]%*%fit$b[(no.basis+extra.knots)*4+1:(no.years-2)],type="l",col=4,xlab="Year",ylab=bquote(kappa[t]),main="Period Effect (Males)")
lines(1961:2016,inv.kappa.trans[,-c(1,2)]%*%fit$b[(no.basis+extra.knots)*4+no.years-2+1:(no.years-2)],type="l",col=1)
lines(1961:2016,inv.kappa.trans[,-c(1,2)]%*%fit.m.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)],type="l",col=4,lty=2)
lines(1961:2016,inv.kappa.trans[,-c(1,2)]%*%fit.m.Scot.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)],type="l",col=1,lty=2)
legend(1967,-0.055,c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,1,4,1),x.intersp=0.2,y.intersp=0.6,seg.len=1)

plot(1857:2015,A2.reduced%*%fit$b[(no.basis+extra.knots)*4+(no.years-2)*2+1:(no.basis-3)],type="l",col=4,xlab="Year of Birth",ylab=bquote(s[gamma](t-x)),main="Cohort Effect (Males)")
lines(1857:2015,c(rep(NA,last.age-last.age.Scot),A2.reduced.Scot%*%fit$b[(no.basis+extra.knots)*4+(no.years-2)*2+no.basis-3+1:ncol(A2.reduced.Scot)]),type="l",col=1)
lines(1857:2015,A2.reduced%*%fit.m.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(no.basis-3)],type="l",col=4,lty=2)
lines(1857:2015,c(rep(NA,last.age-last.age.Scot),A2.reduced.Scot%*%fit.m.Scot.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:ncol(A2.reduced.Scot)]),type="l",col=1,lty=2)
legend(1967,-0.055,c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,1,4,1),x.intersp=0.2,y.intersp=0.6,seg.len=1)


plot(inv.kappa.trans[,-c(1,2)]%*%fit$b[(no.basis+extra.knots)*4+1:(no.years-2)],type="l",col=4)
lines(inv.kappa.trans[,-c(1,2)]%*%fit$b[(no.basis+extra.knots)*4+no.years-2+1:(no.years-2)],type="l",col=2)
lines(inv.kappa.trans[,-c(1,2)]%*%fit.m.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)],type="l",col=4,lty=2)
lines(inv.kappa.trans[,-c(1,2)]%*%fit.m.Scot.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)],type="l",col=2,lty=2)

plot(A2.reduced%*%fit$b[(no.basis+extra.knots)*4+(no.years-2)*2+1:(no.basis-3)],type="l",col=4)
lines(c(rep(0,4),A2.reduced.Scot%*%fit$b[(no.basis+extra.knots)*4+(no.years-2)*2+no.basis-3+1:ncol(A2.reduced.Scot)]),type="l",col=2)
lines(A2.reduced%*%fit.m.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(no.basis-3)],type="l",col=4,lty=2)
lines(c(rep(0,4),A2.reduced.Scot%*%fit.m.Scot.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:ncol(A2.reduced.Scot)]),type="l",col=2,lty=2)

###diff scale
par(mfrow=c(1,2),mar=c(5,6,4,2),cex.lab=1.8,cex.axis=1.3,cex.main=1.5)
plot(1962:2016,diff(inv.kappa.trans[,-c(1,2)]%*%fit$b[(no.basis+extra.knots)*4+1:(no.years-2)]),type="l",col=4,xlab="Year",ylab=as.expression(bquote(kappa[t]^"*")),main="Period Effect on difference scale (Males)")
lines(1962:2016,diff(inv.kappa.trans[,-c(1,2)]%*%fit$b[(no.basis+extra.knots)*4+no.years-2+1:(no.years-2)]),type="l",col=1)
lines(1962:2016,diff(inv.kappa.trans[,-c(1,2)]%*%fit.m.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)]),type="l",col=4,lty=2)
lines(1962:2016,diff(inv.kappa.trans[,-c(1,2)]%*%fit.m.Scot.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)]),type="l",col=1,lty=2)
legend(1970,0.067,c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,1,4,1),y.intersp=0.6,x.intersp=0.2,seg.len=1)

plot(1864:2015,c(rep(NA,last.age-last.age.Scot),diff(A2.reduced.Scot%*%fit$b[(no.basis+extra.knots)*4+(no.years-2)*2+no.basis-3+1:ncol(A2.reduced.Scot)])[-(1:6)]),type="l",col=1,xlab="Year of Birth",ylab=as.expression(bquote(s[gamma]^"*"~"(t-x)")),main="Cohort Effect on difference scale (Males)")
lines(1864:2015,diff(A2.reduced%*%fit$b[(no.basis+extra.knots)*4+(no.years-2)*2+1:(no.basis-3)])[-(1:6)],type="l",col=4)
lines(1864:2015,diff(A2.reduced%*%fit.m.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(no.basis-3)])[-(1:6)],type="l",col=4,lty=2)
lines(1864:2015,c(rep(NA,last.age-last.age.Scot),diff(A2.reduced.Scot%*%fit.m.Scot.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:ncol(A2.reduced.Scot)])[-(1:6)]),type="l",col=1,lty=2)
legend(1840,0.032,c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,1,4,1),y.intersp=0.6,x.intersp=0.2,seg.len=1)

par(mfrow=c(1,2),mar=c(5,6,4,2),cex.lab=1.8,cex.axis=1.3,cex.main=1.5)
plot(1962:2016,diff(inv.kappa.trans[,-c(1,2)]%*%fit.f$b[(no.basis+extra.knots)*4+1:(no.years-2)]),type="l",col=2,xlab="Year",ylab=as.expression(bquote(kappa[t]^"*")),main="Period Effect on difference scale (Females)")
lines(1962:2016,diff(inv.kappa.trans[,-c(1,2)]%*%fit.f$b[(no.basis+extra.knots)*4+no.years-2+1:(no.years-2)]),type="l",col=3)
lines(1962:2016,diff(inv.kappa.trans[,-c(1,2)]%*%fit.f.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)]),type="l",col=2,lty=2)
lines(1962:2016,diff(inv.kappa.trans[,-c(1,2)]%*%fit.f.Scot.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)]),type="l",col=3,lty=2)
legend(1960,-0.06,c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(2,3,2,3),y.intersp=0.6,x.intersp=0.2,seg.len=1)

plot(1864:2015,diff(A2.reduced%*%fit.f$b[(no.basis+extra.knots)*4+(no.years-2)*2+1:(no.basis-3)])[-(1:6)],type="l",col=2,xlab="Year of Birth",ylab=as.expression(bquote(s[gamma]^"*"~"(t-x)")),main="Cohort Effect on difference scale (Females)",ylim=c(-0.015,0.01))
lines(1864:2015,c(rep(NA,last.age-last.age.Scot),diff(A2.reduced.Scot%*%fit.f$b[(no.basis+extra.knots)*4+(no.years-2)*2+no.basis-3+1:ncol(A2.reduced.Scot)])[-(1:6)]),type="l",col=3)
lines(1864:2015,diff(A2.reduced%*%fit.f.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(no.basis-3)])[-(1:6)],type="l",col=2,lty=2)
lines(1864:2015,c(rep(NA,last.age-last.age.Scot),diff(A2.reduced.Scot%*%fit.f.Scot.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:ncol(A2.reduced.Scot)])[-(1:6)]),type="l",col=3,lty=2)
legend(1840,-0.01,c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(2,3,2,3),y.intersp=0.6,x.intersp=0.2,seg.len=1)

######ALL TGT
par(mar=c(5,5,4,2))

plot(A.extra%*%fit.Eng.Scot$b[1:(no.basis+extra.knots)],type="l",col=4,ylim=c(-9,1),cex.lab=1.5,cex.main=1.5,main="Basline Mortality Schedule Extrapolated to 120",xlab="Age",ylab=bquote(s[alpha](x)))
lines(A.extra%*%fit.Eng.Scot$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)],type="l",col=2)
lines(A.extra%*%fit.Eng.Scot$b[(no.basis+extra.knots)*2+1:(no.basis+extra.knots)],type="l",col=1)
lines(A.extra%*%fit.Eng.Scot$b[(no.basis+extra.knots)*3+1:(no.basis+extra.knots)],type="l",col=6)
lines(A.extra%*%fit.m.sep$b[1:(no.basis+extra.knots)],type="l",col=4,lty=2)
lines(A.extra%*%fit.f.sep$b[1:(no.basis+extra.knots)],type="l",col=2,lty=2)
lines(A.extra%*%fit.m.Scot.sep$b[1:(no.basis+extra.knots)],type="l",col=1,lty=2)
lines(A.extra%*%fit.f.Scot.sep$b[1:(no.basis+extra.knots)],type="l",col=6,lty=2)
legend(-5,1,c("EW Males (Joint)","EW Females (Joint)","Scotland Males (Joint)","Scotland Females (Joint)","EW Males (Separate)","EW Females (Separate)","Scotland Males (Separate)","Scotland Females (Separate)"),bty="n",lty=c(1,1,1,1,2,2,2,2),col=c(4,2,1,6,4,2,1,6),x.intersp=0.3,y.intersp=0.7,seg.len=1)

plot(A.extra%*%fit.Eng.Scot$b[(no.basis+extra.knots)*4+1:(no.basis+extra.knots)],type="l",col=4,ylim=c(-2.2,0.3),cex.lab=1.5,cex.main=1.5,main="Age-specific Improvement Rates Extrapolated to 120",xlab="Age",ylab=bquote(s[beta](x)))
lines(A.extra%*%fit.Eng.Scot$b[(no.basis+extra.knots)*5+1:(no.basis+extra.knots)],type="l",col=2)
lines(A.extra%*%fit.Eng.Scot$b[(no.basis+extra.knots)*6+1:(no.basis+extra.knots)],type="l",col=1)
lines(A.extra%*%fit.Eng.Scot$b[(no.basis+extra.knots)*7+1:(no.basis+extra.knots)],type="l",col=3)
lines(A.extra%*%fit.m.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)],type="l",col=4,lty=2)
lines(A.extra%*%fit.f.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)],type="l",col=2,lty=2)
lines(A.extra%*%fit.m.Scot.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)],type="l",col=1,lty=2)
lines(A.extra%*%fit.f.Scot.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)],type="l",col=3,lty=2)
legend(90,-1.3,c("EW Males (Joint)","EW Females (Joint)","Scotland Males (Joint)","Scotland Females (Joint)","EW Males (Separate)","EW Females (Separate)","Scotland Males (Separate)","Scotland Females (Separate)"),bty="n",lty=c(1,1,1,1,2,2,2,2),col=c(4,2,1,6,4,2,1,6),x.intersp=0.3,y.intersp=0.7,seg.len=1)

lines(A.extra%*%fit$b[(no.basis+extra.knots)*2+1:(no.basis+extra.knots)],type="l",col=4,lty=3)
lines(A.extra%*%fit$b[(no.basis+extra.knots)*3+1:(no.basis+extra.knots)],type="l",col=1,lty=3)
lines(A.extra%*%fit.f$b[(no.basis+extra.knots)*2+1:(no.basis+extra.knots)],type="l",col=2,lty=3)
lines(A.extra%*%fit.f$b[(no.basis+extra.knots)*3+1:(no.basis+extra.knots)],type="l",col=6,lty=3)

plot(inv.kappa.trans[,-c(1,2)]%*%fit.f$b[(no.basis+extra.knots)*4+1:(no.years-2)],type="l",col=4)
lines(inv.kappa.trans[,-c(1,2)]%*%fit.f$b[(no.basis+extra.knots)*4+no.years-2+1:(no.years-2)],type="l",col=2)
lines(inv.kappa.trans[,-c(1,2)]%*%fit.f.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)],type="l",col=4,lty=2)
lines(inv.kappa.trans[,-c(1,2)]%*%fit.f.Scot.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)],type="l",col=2,lty=2)

plot(A2.reduced%*%fit.f$b[(no.basis+extra.knots)*4+(no.years-2)*2+1:(no.basis-3)],type="l",col=4)
lines(c(rep(0,4),A2.reduced.Scot%*%fit.f$b[(no.basis+extra.knots)*4+(no.years-2)*2+no.basis-3+1:ncol(A2.reduced.Scot)]),type="l",col=2)
lines(A2.reduced%*%fit.f.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(no.basis-3)],type="l",col=4,lty=2)
lines(c(rep(0,4),A2.reduced.Scot%*%fit.f.Scot.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:ncol(A2.reduced.Scot)]),type="l",col=2,lty=2)

forecast.years=25
A2.extra<-c()
for(j in 1:(no.basis)) {
A2.extra<-cbind(A2.extra,bspline(c(cohort.years,1+diff(cohort.years)[1]*1:100),knots,j))
}
forecast.time<-0.5+diff(timeindex)[1]*forecast.years
forecast.time.50<-0.5+diff(timeindex)[1]*50
forecast.time.100<-0.5+diff(timeindex)[1]*100

####males and females
cohort.years.m<-A2.extra%*%(C.Q.trans%*%fit.Eng$b[(no.basis+extra.knots)*4+(no.years-2)*2+1:(no.basis-3)])
cohort.years.f<-A2.extra%*%(C.Q.trans%*%fit.Eng$b[(no.basis+extra.knots)*4+(no.years-2)*2+no.basis-3+1:(no.basis-3)])
cohort.years.m.Scot<-A2.extra%*%c(0,(C.Q.trans.Scot%*%fit.Scot$b[(no.basis+extra.knots)*4+(no.years-2)*2+1:(ncol(A2.Scot)-3)]))
cohort.years.f.Scot<-A2.extra%*%c(0,(C.Q.trans.Scot%*%fit.Scot$b[(no.basis+extra.knots)*4+(no.years-2)*2+ncol(A2.Scot)-3+1:(ncol(A2.Scot)-3)]))

alpha.m<-A.extra%*%fit.Eng$b[(no.basis+extra.knots)*0+1:(no.basis+extra.knots)]
alpha.f<-A.extra%*%fit.Eng$b[(no.basis+extra.knots)*1+1:(no.basis+extra.knots)]
beta.m<-A.extra%*%fit.Eng$b[(no.basis+extra.knots)*2+1:(no.basis+extra.knots)]
beta.f<-A.extra%*%fit.Eng$b[(no.basis+extra.knots)*3+1:(no.basis+extra.knots)]
kappa.m<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.Eng$b[(no.basis+extra.knots)*4+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
kappa.f<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.Eng$b[(no.basis+extra.knots)*4+no.years-2+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
cohort.m<-cohort.years.m[length(cohort.years)+forecast.years:(forecast.years-120+1)]
cohort.f<-cohort.years.f[length(cohort.years)+forecast.years:(forecast.years-120+1)]

alpha.m.Scot<-A.extra%*%fit.Scot$b[(no.basis+extra.knots)*0+1:(no.basis+extra.knots)]
alpha.f.Scot<-A.extra%*%fit.Scot$b[(no.basis+extra.knots)*1+1:(no.basis+extra.knots)]
beta.m.Scot<-A.extra%*%fit.Scot$b[(no.basis+extra.knots)*2+1:(no.basis+extra.knots)]
beta.f.Scot<-A.extra%*%fit.Scot$b[(no.basis+extra.knots)*3+1:(no.basis+extra.knots)]
kappa.m.Scot<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.Scot$b[(no.basis+extra.knots)*4+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
kappa.f.Scot<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.Scot$b[(no.basis+extra.knots)*4+no.years-2+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
cohort.m.Scot<-cohort.years.m.Scot[length(cohort.years)+forecast.years:(forecast.years-120+1)]
cohort.f.Scot<-cohort.years.f.Scot[length(cohort.years)+forecast.years:(forecast.years-120+1)]

cohort.years.m.sep<-A2.extra%*%(C.Q.trans%*%fit.m.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(no.basis-3)])
cohort.years.m.Scot.sep<-A2.extra%*%c(0,(C.Q.trans.Scot%*%fit.m.Scot.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(ncol(A2.Scot)-3)]))
cohort.years.f.sep<-A2.extra%*%(C.Q.trans%*%fit.f.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(no.basis-3)])
cohort.years.f.Scot.sep<-A2.extra%*%c(0,(C.Q.trans.Scot%*%fit.f.Scot.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(ncol(A2.Scot)-3)]))

alpha.m.sep<-A.extra%*%fit.m.sep$b[1:(no.basis+extra.knots)]
alpha.m.Scot.sep<-A.extra%*%fit.m.Scot.sep$b[1:(no.basis+extra.knots)]
beta.m.sep<-A.extra%*%fit.m.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)]
beta.m.Scot.sep<-A.extra%*%fit.m.Scot.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)]
kappa.m.sep<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.m.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
kappa.m.Scot.sep<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.m.Scot.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
cohort.m.sep<-cohort.years.m.sep[length(cohort.years)+forecast.years:(forecast.years-120+1)]
cohort.m.Scot.sep<-cohort.years.m.Scot.sep[length(cohort.years)+forecast.years:(forecast.years-120+1)]

alpha.f.sep<-A.extra%*%fit.f.sep$b[1:(no.basis+extra.knots)]
alpha.f.Scot.sep<-A.extra%*%fit.f.Scot.sep$b[1:(no.basis+extra.knots)]
beta.f.sep<-A.extra%*%fit.f.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)]
beta.f.Scot.sep<-A.extra%*%fit.f.Scot.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)]
kappa.f.sep<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.f.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
kappa.f.Scot.sep<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.f.Scot.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
cohort.f.sep<-cohort.years.f.sep[length(cohort.years)+forecast.years:(forecast.years-120+1)]
cohort.f.Scot.sep<-cohort.years.f.Scot.sep[length(cohort.years)+forecast.years:(forecast.years-120+1)]

par(mfrow=c(1,3),mar=c(5,5,4,2))
plot(alpha.m+beta.m*forecast.time+kappa.m$pred[forecast.years]+cohort.m,type="l",col=4,cex.lab=1.5,cex.main=1.5,main="25-year ahead EW Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.f+beta.f*forecast.time+kappa.f$pred[forecast.years]+cohort.f,type="l",col=2)
lines(alpha.m.sep+beta.m.sep*forecast.time+kappa.m.sep$pred[forecast.years]+cohort.m.sep,col=4,lty=2)
lines(alpha.f.sep+beta.f.sep*forecast.time+kappa.f.sep$pred[forecast.years]+cohort.f.sep,col=2,lty=2)
legend("topleft",c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,2,4,2))

plot(alpha.m+beta.m*forecast.time.50+kappa.m$pred[50]+cohort.years.m[length(cohort.years)+50:(50-120+1)],type="l",col=4,cex.lab=1.5,cex.main=1.5,main="50-year ahead EW Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.f+beta.f*forecast.time.50+kappa.f$pred[50]+cohort.years.f[length(cohort.years)+50:(50-120+1)],type="l",col=2)
lines(alpha.m.sep+beta.m.sep*forecast.time.50+kappa.m.sep$pred[50]+cohort.years.m.sep[length(cohort.years)+50:(50-120+1)],col=4,lty=2)
lines(alpha.f.sep+beta.f.sep*forecast.time.50+kappa.f.sep$pred[50]+cohort.years.f.sep[length(cohort.years)+50:(50-120+1)],col=2,lty=2)

plot(alpha.m+beta.m*forecast.time.100+kappa.m$pred[100]+cohort.years.m[length(cohort.years)+100:(100-120+1)],type="l",col=4,cex.lab=1.5,cex.main=1.5,main="100-year ahead EW Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.f+beta.f*forecast.time.100+kappa.f$pred[100]+cohort.years.f[length(cohort.years)+100:(100-120+1)],type="l",col=2)
lines(alpha.m.sep+beta.m.sep*forecast.time.100+kappa.m.sep$pred[100]+cohort.years.m.sep[length(cohort.years)+100:(100-120+1)],col=4,lty=2)
lines(alpha.f.sep+beta.f.sep*forecast.time.100+kappa.f.sep$pred[100]+cohort.years.f.sep[length(cohort.years)+100:(100-120+1)],col=2,lty=2)


plot(alpha.m.Scot+beta.m.Scot*forecast.time+kappa.m.Scot$pred[forecast.years]+cohort.m.Scot,type="l",col=4,cex.lab=1.5,cex.main=1.5,main="25-year ahead Scotland Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.f.Scot+beta.f.Scot*forecast.time+kappa.f.Scot$pred[forecast.years]+cohort.f.Scot,col=2)
lines(alpha.m.Scot.sep+beta.m.Scot.sep*forecast.time+kappa.m.Scot.sep$pred[forecast.years]+cohort.m.Scot.sep,col=4,lty=2)
lines(alpha.f.Scot.sep+beta.f.Scot.sep*forecast.time+kappa.f.Scot.sep$pred[forecast.years]+cohort.f.Scot.sep,col=2,lty=2)
legend("topleft",c("Males (Joint)","Females (Joint)","Males (Separate)","Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,2,4,2))

plot(alpha.m.Scot+beta.m.Scot*forecast.time.50+kappa.m.Scot$pred[50]+cohort.years.m.Scot[length(cohort.years)+50:(50-120+1)],type="l",col=4,cex.lab=1.5,cex.main=1.5,main="50-year ahead Scotland Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.f.Scot+beta.f.Scot*forecast.time.50+kappa.f.Scot$pred[50]+cohort.years.f.Scot[length(cohort.years)+50:(50-120+1)],type="l",col=2)
lines(alpha.m.Scot.sep+beta.m.Scot.sep*forecast.time.50+kappa.m.Scot.sep$pred[50]+cohort.years.m.Scot.sep[length(cohort.years)+50:(50-120+1)],col=4,lty=2)
lines(alpha.f.Scot.sep+beta.f.Scot.sep*forecast.time.50+kappa.f.Scot.sep$pred[50]+cohort.years.f.Scot.sep[length(cohort.years)+50:(50-120+1)],col=2,lty=2)

plot(alpha.m.Scot+beta.m.Scot*forecast.time.100+kappa.m.Scot$pred[100]+cohort.years.m.Scot[length(cohort.years)+100:(100-120+1)],type="l",col=4,cex.lab=1.5,cex.main=1.5,main="100-year ahead Scotland Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.f.Scot+beta.f.Scot*forecast.time.100+kappa.f.Scot$pred[100]+cohort.years.f.Scot[length(cohort.years)+100:(100-120+1)],type="l",col=2)
lines(alpha.m.Scot.sep+beta.m.Scot.sep*forecast.time.100+kappa.m.Scot.sep$pred[100]+cohort.years.m.Scot.sep[length(cohort.years)+100:(100-120+1)],col=4,lty=2)
lines(alpha.f.Scot.sep+beta.f.Scot.sep*forecast.time.100+kappa.f.Scot.sep$pred[100]+cohort.years.f.Scot.sep[length(cohort.years)+100:(100-120+1)],col=2,lty=2)

#####Eng and Scot
cohort.years.m<-A2.extra%*%(C.Q.trans%*%fit$b[(no.basis+extra.knots)*4+(no.years-2)*2+1:(no.basis-3)])
cohort.years.m.Scot<-A2.extra%*%c(0,(C.Q.trans.Scot%*%fit$b[(no.basis+extra.knots)*4+(no.years-2)*2+no.basis-3+1:(ncol(A2.Scot)-3)]))
cohort.years.f<-A2.extra%*%(C.Q.trans%*%fit.f$b[(no.basis+extra.knots)*4+(no.years-2)*2+1:(no.basis-3)])
cohort.years.f.Scot<-A2.extra%*%c(0,(C.Q.trans.Scot%*%fit.f$b[(no.basis+extra.knots)*4+(no.years-2)*2+no.basis-3+1:(ncol(A2.Scot)-3)]))

alpha.m<-A.extra%*%fit$b[(no.basis+extra.knots)*0+1:(no.basis+extra.knots)]
alpha.m.Scot<-A.extra%*%fit$b[(no.basis+extra.knots)*1+1:(no.basis+extra.knots)]
beta.m<-A.extra%*%fit$b[(no.basis+extra.knots)*2+1:(no.basis+extra.knots)]
beta.m.Scot<-A.extra%*%fit$b[(no.basis+extra.knots)*3+1:(no.basis+extra.knots)]
kappa.m<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit$b[(no.basis+extra.knots)*4+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
kappa.m.Scot<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit$b[(no.basis+extra.knots)*4+no.years-2+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
cohort.m<-cohort.years.m[length(cohort.years)+forecast.years:(forecast.years-120+1)]
cohort.m.Scot<-cohort.years.m.Scot[length(cohort.years)+forecast.years:(forecast.years-120+1)]

alpha.f<-A.extra%*%fit.f$b[(no.basis+extra.knots)*0+1:(no.basis+extra.knots)]
alpha.f.Scot<-A.extra%*%fit.f$b[(no.basis+extra.knots)*1+1:(no.basis+extra.knots)]
beta.f<-A.extra%*%fit.f$b[(no.basis+extra.knots)*2+1:(no.basis+extra.knots)]
beta.f.Scot<-A.extra%*%fit.f$b[(no.basis+extra.knots)*3+1:(no.basis+extra.knots)]
kappa.f<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.f$b[(no.basis+extra.knots)*4+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
kappa.f.Scot<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.f$b[(no.basis+extra.knots)*4+no.years-2+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
cohort.f<-cohort.years.f[length(cohort.years)+forecast.years:(forecast.years-120+1)]
cohort.f.Scot<-cohort.years.f.Scot[length(cohort.years)+forecast.years:(forecast.years-120+1)]

cohort.years.m.sep<-A2.extra%*%(C.Q.trans%*%fit.m.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(no.basis-3)])
cohort.years.m.Scot.sep<-A2.extra%*%c(0,(C.Q.trans.Scot%*%fit.m.Scot.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(ncol(A2.Scot)-3)]))
cohort.years.f.sep<-A2.extra%*%(C.Q.trans%*%fit.f.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(no.basis-3)])
cohort.years.f.Scot.sep<-A2.extra%*%c(0,(C.Q.trans.Scot%*%fit.f.Scot.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(ncol(A2.Scot)-3)]))

alpha.m.sep<-A.extra%*%fit.m.sep$b[1:(no.basis+extra.knots)]
alpha.m.Scot.sep<-A.extra%*%fit.m.Scot.sep$b[1:(no.basis+extra.knots)]
beta.m.sep<-A.extra%*%fit.m.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)]
beta.m.Scot.sep<-A.extra%*%fit.m.Scot.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)]
kappa.m.sep<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.m.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
kappa.m.Scot.sep<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.m.Scot.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
cohort.m.sep<-cohort.years.m.sep[length(cohort.years)+forecast.years:(forecast.years-120+1)]
cohort.m.Scot.sep<-cohort.years.m.Scot.sep[length(cohort.years)+forecast.years:(forecast.years-120+1)]

alpha.f.sep<-A.extra%*%fit.f.sep$b[1:(no.basis+extra.knots)]
alpha.f.Scot.sep<-A.extra%*%fit.f.Scot.sep$b[1:(no.basis+extra.knots)]
beta.f.sep<-A.extra%*%fit.f.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)]
beta.f.Scot.sep<-A.extra%*%fit.f.Scot.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)]
kappa.f.sep<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.f.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
kappa.f.Scot.sep<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.f.Scot.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
cohort.f.sep<-cohort.years.f.sep[length(cohort.years)+forecast.years:(forecast.years-120+1)]
cohort.f.Scot.sep<-cohort.years.f.Scot.sep[length(cohort.years)+forecast.years:(forecast.years-120+1)]



par(mfrow=c(1,3),mar=c(5,5,4,2))
plot(alpha.m+beta.m*forecast.time+kappa.m$pred[forecast.years]+cohort.m,type="l",col=4,cex.lab=1.5,cex.main=1.5,main="25-year ahead Male Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.m.Scot+beta.m.Scot*forecast.time+kappa.m.Scot$pred[forecast.years]+cohort.m.Scot,type="l",col=1)
lines(alpha.m.sep+beta.m.sep*forecast.time+kappa.m.sep$pred[forecast.years]+cohort.m.sep,col=4,lty=2)
lines(alpha.m.Scot.sep+beta.m.Scot.sep*forecast.time+kappa.m.Scot.sep$pred[forecast.years]+cohort.m.Scot.sep,col=1,lty=2)
legend("topleft",c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,1,4,1))

plot(alpha.m+beta.m*forecast.time.50+kappa.m$pred[50]+cohort.years.m[length(cohort.years)+50:(50-120+1)],type="l",col=4,cex.lab=1.5,cex.main=1.5,main="50-year ahead Male Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.m.Scot+beta.m.Scot*forecast.time.50+kappa.m.Scot$pred[50]+cohort.years.m.Scot[length(cohort.years)+50:(50-120+1)],type="l",col=1)
lines(alpha.m.sep+beta.m.sep*forecast.time.50+kappa.m.sep$pred[50]+cohort.years.m.sep[length(cohort.years)+50:(50-120+1)],col=4,lty=2)
lines(alpha.m.Scot.sep+beta.m.Scot.sep*forecast.time.50+kappa.m.Scot.sep$pred[50]+cohort.years.m.Scot.sep[length(cohort.years)+50:(50-120+1)],col=1,lty=2)

plot(alpha.m+beta.m*forecast.time.100+kappa.m$pred[100]+cohort.years.m[length(cohort.years)+100:(100-120+1)],type="l",col=4,cex.lab=1.5,cex.main=1.5,main="100-year ahead Male Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.m.Scot+beta.m.Scot*forecast.time.100+kappa.m.Scot$pred[100]+cohort.years.m.Scot[length(cohort.years)+100:(100-120+1)],type="l",col=1)
lines(alpha.m.sep+beta.m.sep*forecast.time.100+kappa.m.sep$pred[100]+cohort.years.m.sep[length(cohort.years)+100:(100-120+1)],col=4,lty=2)
lines(alpha.m.Scot.sep+beta.m.Scot.sep*forecast.time.100+kappa.m.Scot.sep$pred[100]+cohort.years.m.Scot.sep[length(cohort.years)+100:(100-120+1)],col=1,lty=2)


plot(alpha.f+beta.f*forecast.time+kappa.f$pred[forecast.years]+cohort.f,type="l",col=2,cex.lab=1.5,cex.main=1.5,main="25-year ahead Female Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.f.Scot+beta.f.Scot*forecast.time+kappa.f.Scot$pred[forecast.years]+cohort.f.Scot,col=6)
lines(alpha.f.sep+beta.f.sep*forecast.time+kappa.f.sep$pred[forecast.years]+cohort.f.sep,col=2,lty=2)
lines(alpha.f.Scot.sep+beta.f.Scot.sep*forecast.time+kappa.f.Scot.sep$pred[forecast.years]+cohort.f.Scot.sep,col=6,lty=2)
legend("topleft",c("EW (Joint)","Scotland (Joint)","EW (Separate)","Scotland (Separate)"),bty="n",lty=c(1,1,2,2),col=c(2,6,2,6))

plot(alpha.f+beta.f*forecast.time.50+kappa.f$pred[50]+cohort.years.f[length(cohort.years)+50:(50-120+1)],type="l",col=2,cex.lab=1.5,cex.main=1.5,main="50-year ahead Female Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.f.Scot+beta.f.Scot*forecast.time.50+kappa.f.Scot$pred[50]+cohort.years.f.Scot[length(cohort.years)+50:(50-120+1)],type="l",col=6)
lines(alpha.f.sep+beta.f.sep*forecast.time.50+kappa.f.sep$pred[50]+cohort.years.f.sep[length(cohort.years)+50:(50-120+1)],col=2,lty=2)
lines(alpha.f.Scot.sep+beta.f.Scot.sep*forecast.time.50+kappa.f.Scot.sep$pred[50]+cohort.years.f.Scot.sep[length(cohort.years)+50:(50-120+1)],col=6,lty=2)

plot(alpha.f+beta.f*forecast.time.100+kappa.f$pred[100]+cohort.years.f[length(cohort.years)+100:(100-120+1)],type="l",col=2,cex.lab=1.5,cex.main=1.5,main="100-year ahead Female Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.f.Scot+beta.f.Scot*forecast.time.100+kappa.f.Scot$pred[100]+cohort.years.f.Scot[length(cohort.years)+100:(100-120+1)],type="l",col=6)
lines(alpha.f.sep+beta.f.sep*forecast.time.100+kappa.f.sep$pred[100]+cohort.years.f.sep[length(cohort.years)+100:(100-120+1)],col=2,lty=2)
lines(alpha.f.Scot.sep+beta.f.Scot.sep*forecast.time.100+kappa.f.Scot.sep$pred[100]+cohort.years.f.Scot.sep[length(cohort.years)+100:(100-120+1)],col=6,lty=2)


####ALL TGT
cohort.years.m<-A2.extra%*%(C.Q.trans%*%fit.Eng.Scot$b[(no.basis+extra.knots)*8+(no.years-2)*4+1:(no.basis-3)])
cohort.years.m.Scot<-A2.extra%*%c(0,(C.Q.trans.Scot%*%fit.Eng.Scot$b[(no.basis+extra.knots)*8+(no.years-2)*4+(no.basis-3)*2+1:(ncol(A2.Scot)-3)]))
cohort.years.f<-A2.extra%*%(C.Q.trans%*%fit.Eng.Scot$b[(no.basis+extra.knots)*8+(no.years-2)*4+no.basis-3+1:(no.basis-3)])
cohort.years.f.Scot<-A2.extra%*%c(0,(C.Q.trans.Scot%*%fit.Eng.Scot$b[(no.basis+extra.knots)*8+(no.years-2)*4+(no.basis-3)*2+ncol(A2.Scot)-3+1:(ncol(A2.Scot)-3)]))

alpha.m<-A.extra%*%fit.Eng.Scot$b[(no.basis+extra.knots)*0+1:(no.basis+extra.knots)]
alpha.m.Scot<-A.extra%*%fit.Eng.Scot$b[(no.basis+extra.knots)*2+1:(no.basis+extra.knots)]
beta.m<-A.extra%*%fit.Eng.Scot$b[(no.basis+extra.knots)*4+1:(no.basis+extra.knots)]
beta.m.Scot<-A.extra%*%fit.Eng.Scot$b[(no.basis+extra.knots)*6+1:(no.basis+extra.knots)]
kappa.m<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.Eng.Scot$b[(no.basis+extra.knots)*8+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
kappa.m.Scot<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.Eng.Scot$b[(no.basis+extra.knots)*8+(no.years-2)*2+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
cohort.m<-cohort.years.m[length(cohort.years)+forecast.years:(forecast.years-120+1)]
cohort.m.Scot<-cohort.years.m.Scot[length(cohort.years)+forecast.years:(forecast.years-120+1)]

alpha.f<-A.extra%*%fit.Eng.Scot$b[(no.basis+extra.knots)*1+1:(no.basis+extra.knots)]
alpha.f.Scot<-A.extra%*%fit.Eng.Scot$b[(no.basis+extra.knots)*3+1:(no.basis+extra.knots)]
beta.f<-A.extra%*%fit.Eng.Scot$b[(no.basis+extra.knots)*5+1:(no.basis+extra.knots)]
beta.f.Scot<-A.extra%*%fit.Eng.Scot$b[(no.basis+extra.knots)*7+1:(no.basis+extra.knots)]
kappa.f<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.Eng.Scot$b[(no.basis+extra.knots)*8+no.years-2+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
kappa.f.Scot<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.Eng.Scot$b[(no.basis+extra.knots)*8+(no.years-2)*3+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
cohort.f<-cohort.years.f[length(cohort.years)+forecast.years:(forecast.years-120+1)]
cohort.f.Scot<-cohort.years.f.Scot[length(cohort.years)+forecast.years:(forecast.years-120+1)]

cohort.years.m.sep<-A2.extra%*%(C.Q.trans%*%fit.m.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(no.basis-3)])
cohort.years.m.Scot.sep<-A2.extra%*%c(0,(C.Q.trans.Scot%*%fit.m.Scot.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(ncol(A2.Scot)-3)]))
cohort.years.f.sep<-A2.extra%*%(C.Q.trans%*%fit.f.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(no.basis-3)])
cohort.years.f.Scot.sep<-A2.extra%*%c(0,(C.Q.trans.Scot%*%fit.f.Scot.sep$b[(no.basis+extra.knots)*2+(no.years-2)+1:(ncol(A2.Scot)-3)]))

alpha.m.sep<-A.extra%*%fit.m.sep$b[1:(no.basis+extra.knots)]
alpha.m.Scot.sep<-A.extra%*%fit.m.Scot.sep$b[1:(no.basis+extra.knots)]
beta.m.sep<-A.extra%*%fit.m.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)]
beta.m.Scot.sep<-A.extra%*%fit.m.Scot.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)]
kappa.m.sep<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.m.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
kappa.m.Scot.sep<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.m.Scot.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
cohort.m.sep<-cohort.years.m.sep[length(cohort.years)+forecast.years:(forecast.years-120+1)]
cohort.m.Scot.sep<-cohort.years.m.Scot.sep[length(cohort.years)+forecast.years:(forecast.years-120+1)]

alpha.f.sep<-A.extra%*%fit.f.sep$b[1:(no.basis+extra.knots)]
alpha.f.Scot.sep<-A.extra%*%fit.f.Scot.sep$b[1:(no.basis+extra.knots)]
beta.f.sep<-A.extra%*%fit.f.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)]
beta.f.Scot.sep<-A.extra%*%fit.f.Scot.sep$b[(no.basis+extra.knots)+1:(no.basis+extra.knots)]
kappa.f.sep<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.f.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
kappa.f.Scot.sep<-predict(arima(inv.kappa.trans[,-c(1,2)]%*%fit.f.Scot.sep$b[(no.basis+extra.knots)*2+1:(no.years-2)],order=c(1,0,0),include.mean=F),n.ahead=100)
cohort.f.sep<-cohort.years.f.sep[length(cohort.years)+forecast.years:(forecast.years-120+1)]
cohort.f.Scot.sep<-cohort.years.f.Scot.sep[length(cohort.years)+forecast.years:(forecast.years-120+1)]

par(mfrow=c(1,3),mar=c(5,5,4,2))
plot(alpha.m+beta.m*forecast.time+kappa.m$pred[forecast.years]+cohort.m,type="l",col=4,cex.lab=1.5,cex.main=1.5,main="25-year ahead Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.f+beta.f*forecast.time+kappa.f$pred[forecast.years]+cohort.f,type="l",col=2,main="Females")
lines(alpha.m.Scot+beta.m.Scot*forecast.time+kappa.m.Scot$pred[forecast.years]+cohort.m.Scot,col=1)
lines(alpha.f.Scot+beta.f.Scot*forecast.time+kappa.f.Scot$pred[forecast.years]+cohort.f.Scot,col=6)
lines(alpha.m.sep+beta.m.sep*forecast.time+kappa.m.sep$pred[forecast.years]+cohort.m.sep,col=4,lty=2)
lines(alpha.f.sep+beta.f.sep*forecast.time+kappa.f.sep$pred[forecast.years]+cohort.f.sep,col=2,lty=2)
lines(alpha.m.Scot.sep+beta.m.Scot.sep*forecast.time+kappa.m.Scot.sep$pred[forecast.years]+cohort.m.Scot.sep,col=1,lty=2)
lines(alpha.f.Scot.sep+beta.f.Scot.sep*forecast.time+kappa.f.Scot.sep$pred[forecast.years]+cohort.f.Scot.sep,col=6,lty=2)
legend(-15,1,c("EW Males (Joint)","EW Females (Joint)", "Scotland Males (Joint)","Scotland Females (Joint)","EW Males (Separate)","EW Females (Separate)", "Scotland Males (Separate)","Scotland Females (Separate)"),bty="n",lty=c(1,1,1,1,2,2,2,2),col=c(4,2,1,6,4,2,1,6),x.intersp=0.3)

plot(alpha.m+beta.m*forecast.time.50+kappa.m$pred[50]+cohort.years.m[length(cohort.years)+50:(50-120+1)],type="l",col=4,cex.lab=1.5,cex.main=1.5,main="50-year ahead Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.f+beta.f*forecast.time.50+kappa.f$pred[50]+cohort.years.f[length(cohort.years)+50:(50-120+1)],type="l",col=2,main="Females")
lines(alpha.m.Scot+beta.m.Scot*forecast.time.50+kappa.m.Scot$pred[50]+cohort.years.m.Scot[length(cohort.years)+50:(50-120+1)],col=1)
lines(alpha.f.Scot+beta.f.Scot*forecast.time.50+kappa.f.Scot$pred[50]+cohort.years.f.Scot[length(cohort.years)+50:(50-120+1)],col=6)
lines(alpha.m.sep+beta.m.sep*forecast.time.50+kappa.m.sep$pred[50]+cohort.years.m.sep[length(cohort.years)+50:(50-120+1)],col=4,lty=2)
lines(alpha.f.sep+beta.f.sep*forecast.time.50+kappa.f.sep$pred[50]+cohort.years.f.sep[length(cohort.years)+50:(50-120+1)],col=2,lty=2)
lines(alpha.m.Scot.sep+beta.m.Scot.sep*forecast.time.50+kappa.m.Scot.sep$pred[50]+cohort.years.m.Scot.sep[length(cohort.years)+50:(50-120+1)],col=1,lty=2)
lines(alpha.f.Scot.sep+beta.f.Scot.sep*forecast.time.50+kappa.f.Scot.sep$pred[50]+cohort.years.f.Scot.sep[length(cohort.years)+50:(50-120+1)],col=6,lty=2)

plot(alpha.m+beta.m*forecast.time.100+kappa.m$pred[100]+cohort.years.m[length(cohort.years)+100:(100-120+1)],type="l",col=4,cex.lab=1.5,cex.main=1.5,main="100-year ahead Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.f+beta.f*forecast.time.100+kappa.f$pred[100]+cohort.years.f[length(cohort.years)+100:(100-120+1)],type="l",col=2,main="Females")
lines(alpha.m.Scot+beta.m.Scot*forecast.time.100+kappa.m.Scot$pred[100]+cohort.years.m.Scot[length(cohort.years)+100:(100-120+1)],col=1)
lines(alpha.f.Scot+beta.f.Scot*forecast.time.100+kappa.f.Scot$pred[100]+cohort.years.f.Scot[length(cohort.years)+100:(100-120+1)],col=6)
lines(alpha.m.sep+beta.m.sep*forecast.time.100+kappa.m.sep$pred[100]+cohort.years.m.sep[length(cohort.years)+100:(100-120+1)],col=4,lty=2)
lines(alpha.f.sep+beta.f.sep*forecast.time.100+kappa.f.sep$pred[100]+cohort.years.f.sep[length(cohort.years)+100:(100-120+1)],col=2,lty=2)
lines(alpha.m.Scot.sep+beta.m.Scot.sep*forecast.time.100+kappa.m.Scot.sep$pred[100]+cohort.years.m.Scot.sep[length(cohort.years)+100:(100-120+1)],col=1,lty=2)
lines(alpha.f.Scot.sep+beta.f.Scot.sep*forecast.time.100+kappa.f.Scot.sep$pred[100]+cohort.years.f.Scot.sep[length(cohort.years)+100:(100-120+1)],col=6,lty=2)

####all tgt just males and females
###EW
par(mfrow=c(1,3),mar=c(5,5,4,2))
plot(alpha.m+beta.m*forecast.time+kappa.m$pred[forecast.years]+cohort.m,type="l",col=4,cex.lab=1.5,cex.main=1.5,main="25-year ahead EW Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.f+beta.f*forecast.time+kappa.f$pred[forecast.years]+cohort.f,type="l",col=2,main="Females")
lines(alpha.m.sep+beta.m.sep*forecast.time+kappa.m.sep$pred[forecast.years]+cohort.m.sep,col=4,lty=2)
lines(alpha.f.sep+beta.f.sep*forecast.time+kappa.f.sep$pred[forecast.years]+cohort.f.sep,col=2,lty=2)
legend(-15,1,c("EW Males (Joint)","EW Females (Joint)","EW Males (Separate)","EW Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,2,4,2),x.intersp=0.3)

plot(alpha.m+beta.m*forecast.time.50+kappa.m$pred[50]+cohort.years.m[length(cohort.years)+50:(50-120+1)],type="l",col=4,cex.lab=1.5,cex.main=1.5,main="50-year ahead EW Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.f+beta.f*forecast.time.50+kappa.f$pred[50]+cohort.years.f[length(cohort.years)+50:(50-120+1)],type="l",col=2,main="Females")
lines(alpha.m.sep+beta.m.sep*forecast.time.50+kappa.m.sep$pred[50]+cohort.years.m.sep[length(cohort.years)+50:(50-120+1)],col=4,lty=2)
lines(alpha.f.sep+beta.f.sep*forecast.time.50+kappa.f.sep$pred[50]+cohort.years.f.sep[length(cohort.years)+50:(50-120+1)],col=2,lty=2)

plot(alpha.m+beta.m*forecast.time.100+kappa.m$pred[100]+cohort.years.m[length(cohort.years)+100:(100-120+1)],type="l",col=4,cex.lab=1.5,cex.main=1.5,main="100-year ahead EW Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.f+beta.f*forecast.time.100+kappa.f$pred[100]+cohort.years.f[length(cohort.years)+100:(100-120+1)],type="l",col=2,main="Females")
lines(alpha.m.sep+beta.m.sep*forecast.time.100+kappa.m.sep$pred[100]+cohort.years.m.sep[length(cohort.years)+100:(100-120+1)],col=4,lty=2)
lines(alpha.f.sep+beta.f.sep*forecast.time.100+kappa.f.sep$pred[100]+cohort.years.f.sep[length(cohort.years)+100:(100-120+1)],col=2,lty=2)

###Scotland
par(mfrow=c(1,3),mar=c(5,5,4,2))
plot(alpha.m.Scot+beta.m.Scot*forecast.time+kappa.m.Scot$pred[forecast.years]+cohort.m.Scot,col=1,cex.lab=1.5,cex.main=1.5,main="25-year ahead Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])),type="l")
lines(alpha.f.Scot+beta.f.Scot*forecast.time+kappa.f.Scot$pred[forecast.years]+cohort.f.Scot,col=6)
lines(alpha.m.Scot.sep+beta.m.Scot.sep*forecast.time+kappa.m.Scot.sep$pred[forecast.years]+cohort.m.Scot.sep,col=1,lty=2)
lines(alpha.f.Scot.sep+beta.f.Scot.sep*forecast.time+kappa.f.Scot.sep$pred[forecast.years]+cohort.f.Scot.sep,col=6,lty=2)
legend(-15,1,c("Scotland Males (Joint)","Scotland Females (Joint)", "Scotland Males (Separate)","Scotland Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(1,6,1,6),x.intersp=0.3)

plot(alpha.m.Scot+beta.m.Scot*forecast.time.50+kappa.m.Scot$pred[50]+cohort.years.m.Scot[length(cohort.years)+50:(50-120+1)],col=1,cex.lab=1.5,cex.main=1.5,main="50-year ahead Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])),type="l")
lines(alpha.f.Scot+beta.f.Scot*forecast.time.50+kappa.f.Scot$pred[50]+cohort.years.f.Scot[length(cohort.years)+50:(50-120+1)],col=6)
lines(alpha.m.Scot.sep+beta.m.Scot.sep*forecast.time.50+kappa.m.Scot.sep$pred[50]+cohort.years.m.Scot.sep[length(cohort.years)+50:(50-120+1)],col=1,lty=2)
lines(alpha.f.Scot.sep+beta.f.Scot.sep*forecast.time.50+kappa.f.Scot.sep$pred[50]+cohort.years.f.Scot.sep[length(cohort.years)+50:(50-120+1)],col=6,lty=2)

plot(alpha.m.Scot+beta.m.Scot*forecast.time.100+kappa.m.Scot$pred[100]+cohort.years.m.Scot[length(cohort.years)+100:(100-120+1)],col=1,cex.lab=1.5,cex.main=1.5,main="100-year ahead Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])),type="l")
lines(alpha.f.Scot+beta.f.Scot*forecast.time.100+kappa.f.Scot$pred[100]+cohort.years.f.Scot[length(cohort.years)+100:(100-120+1)],col=6)
lines(alpha.m.Scot.sep+beta.m.Scot.sep*forecast.time.100+kappa.m.Scot.sep$pred[100]+cohort.years.m.Scot.sep[length(cohort.years)+100:(100-120+1)],col=1,lty=2)
lines(alpha.f.Scot.sep+beta.f.Scot.sep*forecast.time.100+kappa.f.Scot.sep$pred[100]+cohort.years.f.Scot.sep[length(cohort.years)+100:(100-120+1)],col=6,lty=2)


####all tgt just EW and Scot
	####males
par(mfrow=c(1,3),mar=c(5,5,4,2))
plot(alpha.m+beta.m*forecast.time+kappa.m$pred[forecast.years]+cohort.m,type="l",col=4,cex.lab=1.5,cex.main=1.5,main="25-year ahead Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.m.Scot+beta.m.Scot*forecast.time+kappa.m.Scot$pred[forecast.years]+cohort.m.Scot,col=1)
lines(alpha.m.sep+beta.m.sep*forecast.time+kappa.m.sep$pred[forecast.years]+cohort.m.sep,col=4,lty=2)
lines(alpha.m.Scot.sep+beta.m.Scot.sep*forecast.time+kappa.m.Scot.sep$pred[forecast.years]+cohort.m.Scot.sep,col=1,lty=2)
legend(-15,1,c("EW Males (Joint)", "Scotland Males (Joint)","EW Males (Separate)","Scotland Males (Separate)"),bty="n",lty=c(1,1,2,2),col=c(4,1,4,1),x.intersp=0.3)

plot(alpha.m+beta.m*forecast.time.50+kappa.m$pred[50]+cohort.years.m[length(cohort.years)+50:(50-120+1)],type="l",col=4,cex.lab=1.5,cex.main=1.5,main="50-year ahead Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.m.Scot+beta.m.Scot*forecast.time.50+kappa.m.Scot$pred[50]+cohort.years.m.Scot[length(cohort.years)+50:(50-120+1)],col=1)
lines(alpha.m.sep+beta.m.sep*forecast.time.50+kappa.m.sep$pred[50]+cohort.years.m.sep[length(cohort.years)+50:(50-120+1)],col=4,lty=2)
lines(alpha.m.Scot.sep+beta.m.Scot.sep*forecast.time.50+kappa.m.Scot.sep$pred[50]+cohort.years.m.Scot.sep[length(cohort.years)+50:(50-120+1)],col=1,lty=2)

plot(alpha.m+beta.m*forecast.time.100+kappa.m$pred[100]+cohort.years.m[length(cohort.years)+100:(100-120+1)],type="l",col=4,cex.lab=1.5,cex.main=1.5,main="100-year ahead Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.m.Scot+beta.m.Scot*forecast.time.100+kappa.m.Scot$pred[100]+cohort.years.m.Scot[length(cohort.years)+100:(100-120+1)],col=1)
lines(alpha.m.sep+beta.m.sep*forecast.time.100+kappa.m.sep$pred[100]+cohort.years.m.sep[length(cohort.years)+100:(100-120+1)],col=4,lty=2)
lines(alpha.m.Scot.sep+beta.m.Scot.sep*forecast.time.100+kappa.m.Scot.sep$pred[100]+cohort.years.m.Scot.sep[length(cohort.years)+100:(100-120+1)],col=1,lty=2)


	###females
par(mfrow=c(1,3),mar=c(5,5,4,2))
plot(alpha.f+beta.f*forecast.time+kappa.f$pred[forecast.years]+cohort.f,type="l",col=2,cex.lab=1.5,cex.main=1.5,main="25-year ahead Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.f.Scot+beta.f.Scot*forecast.time+kappa.f.Scot$pred[forecast.years]+cohort.f.Scot,col=6)
lines(alpha.f.sep+beta.f.sep*forecast.time+kappa.f.sep$pred[forecast.years]+cohort.f.sep,col=2,lty=2)
lines(alpha.f.Scot.sep+beta.f.Scot.sep*forecast.time+kappa.f.Scot.sep$pred[forecast.years]+cohort.f.Scot.sep,col=6,lty=2)
legend(-15,1,c("EW Females (Joint)","Scotland Females (Joint)","EW Females (Separate)","Scotland Females (Separate)"),bty="n",lty=c(1,1,2,2),col=c(2,6,2,6),x.intersp=0.3)

plot(alpha.f+beta.f*forecast.time.50+kappa.f$pred[50]+cohort.years.f[length(cohort.years)+50:(50-120+1)],type="l",col=2,cex.lab=1.5,cex.main=1.5,main="50-year ahead Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.f.Scot+beta.f.Scot*forecast.time.50+kappa.f.Scot$pred[50]+cohort.years.f.Scot[length(cohort.years)+50:(50-120+1)],col=6)
lines(alpha.f.sep+beta.f.sep*forecast.time.50+kappa.f.sep$pred[50]+cohort.years.f.sep[length(cohort.years)+50:(50-120+1)],col=2,lty=2)
lines(alpha.f.Scot.sep+beta.f.Scot.sep*forecast.time.50+kappa.f.Scot.sep$pred[50]+cohort.years.f.Scot.sep[length(cohort.years)+50:(50-120+1)],col=6,lty=2)

plot(alpha.f+beta.f*forecast.time.100+kappa.f$pred[100]+cohort.years.f[length(cohort.years)+100:(100-120+1)],type="l",col=2,cex.lab=1.5,cex.main=1.5,main="100-year ahead Mortality Forecasts \n(on Log Scale)",xlab="Age",ylab=bquote(log(m[x])))
lines(alpha.f.Scot+beta.f.Scot*forecast.time.100+kappa.f.Scot$pred[100]+cohort.years.f.Scot[length(cohort.years)+100:(100-120+1)],col=6)
lines(alpha.f.sep+beta.f.sep*forecast.time.100+kappa.f.sep$pred[100]+cohort.years.f.sep[length(cohort.years)+100:(100-120+1)],col=2,lty=2)
lines(alpha.f.Scot.sep+beta.f.Scot.sep*forecast.time.100+kappa.f.Scot.sep$pred[100]+cohort.years.f.Scot.sep[length(cohort.years)+100:(100-120+1)],col=6,lty=2)


####Graduation
no.basis=40;extra.knots=0;last.age=104
knots<-seq(0,1,length=no.basis-2)
dk<-knots[2]-knots[1]	
knots<-c(knots[1]-dk*(3:1),knots,knots[no.basis-2]+dk*(1:(3+extra.knots)))
age<-seq(0,1,length=last.age)

A<-c()
for(j in 1:no.basis) {
A<-cbind(A,bspline(age,knots,j))
}

P<-diff(diag(no.basis+extra.knots),differences=2)
PD<-t(c(1,-1))%x%diag(no.basis+extra.knots)

m.1980<-newton.raphson.bic(data=dm[2:(last.age+1),1980-1961+1],offset=log(Em[2:(last.age+1),1980-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(0),lambda=c(4),typsize=c(1),fscale=1,stepmax=10)
m.2010<-newton.raphson.bic(data=dm[2:(last.age+1),2010-1961+1],offset=log(Em[2:(last.age+1),2010-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(0),lambda=c(4),typsize=c(1),fscale=1,stepmax=10)
m.2011<-newton.raphson.bic(data=dm[2:(last.age+1),2011-1961+1],offset=log(Em[2:(last.age+1),2011-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(0),lambda=c(4),typsize=c(1),fscale=1,stepmax=10)
m.2012<-newton.raphson.bic(data=dm[2:(last.age+1),2012-1961+1],offset=log(Em[2:(last.age+1),2012-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(0),lambda=c(4),typsize=c(1),fscale=1,stepmax=10)
m.2017<-newton.raphson.bic(data=dm[2:(last.age+1),2017-1961+1],offset=log(Em[2:(last.age+1),2017-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(0),lambda=c(4),typsize=c(1),fscale=1,stepmax=10)
m.1980.expo<-newton.raphson.bic(data=dm[2:(last.age+1),1980-1961+1],offset=log(Em[2:(last.age+1),1980-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(1),lambda=c(2,8),typsize=c(1,1),fscale=1,stepmax=10)
m.2007.expo<-newton.raphson.bic(data=dm[2:(last.age+1),2007-1961+1],offset=log(Em[2:(last.age+1),2007-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(1),lambda=c(2,8),typsize=c(1,1),fscale=1,stepmax=5)
m.2010.expo<-newton.raphson.bic(data=dm[2:(last.age+1),2010-1961+1],offset=log(Em[2:(last.age+1),2010-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(1),lambda=c(2,8),typsize=c(1,1),fscale=1,stepmax=5)
m.2011.expo<-newton.raphson.bic(data=dm[2:(last.age+1),2011-1961+1],offset=log(Em[2:(last.age+1),2011-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(1),lambda=c(2,8),typsize=c(1,1),fscale=1,stepmax=5)
m.2012.expo<-newton.raphson.bic(data=dm[2:(last.age+1),2012-1961+1],offset=log(Em[2:(last.age+1),2012-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(1),lambda=c(2,8),typsize=c(1,1),fscale=1,stepmax=5)
m.2017.expo<-newton.raphson.bic(data=dm[2:(last.age+1),2017-1961+1],offset=log(Em[2:(last.age+1),2017-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(1),lambda=c(2,8),typsize=c(1,1),fscale=1,stepmax=5)
f.2007.expo<-newton.raphson.bic(data=df[2:(last.age+1),2007-1961+1],offset=log(Ef[2:(last.age+1),2007-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(1),lambda=c(2,8),typsize=c(1,1),fscale=1,stepmax=5)
f.2010<-newton.raphson.bic(data=df[2:(last.age+1),2010-1961+1],offset=log(Ef[2:(last.age+1),2010-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(0),lambda=c(2),typsize=c(1),fscale=1,stepmax=2)
f.2011<-newton.raphson.bic(data=df[2:(last.age+1),2011-1961+1],offset=log(Ef[2:(last.age+1),2011-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(0),lambda=c(4),typsize=c(1),fscale=1,stepmax=10)
f.2012<-newton.raphson.bic(data=df[2:(last.age+1),2012-1961+1],offset=log(Ef[2:(last.age+1),2012-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(0),lambda=c(4),typsize=c(1),fscale=1,stepmax=10)
f.1980.expo<-newton.raphson.bic(data=df[2:(last.age+1),1980-1961+1],offset=log(Ef[2:(last.age+1),1980-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(1),lambda=c(2,8),typsize=c(1,1),fscale=1,stepmax=5)
f.2010.expo<-newton.raphson.bic(data=df[2:(last.age+1),2010-1961+1],offset=log(Ef[2:(last.age+1),2010-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(1),lambda=c(2,8),typsize=c(1,1),fscale=1,stepmax=5)
f.2011.expo<-newton.raphson.bic(data=df[2:(last.age+1),2011-1961+1],offset=log(Ef[2:(last.age+1),2011-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(1),lambda=c(2,8),typsize=c(1,1),fscale=1,stepmax=5)
f.2012.expo<-newton.raphson.bic(data=df[2:(last.age+1),2012-1961+1],offset=log(Ef[2:(last.age+1),2012-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(1),lambda=c(2,8),typsize=c(1,1),fscale=1,stepmax=5)
f.2017.expo<-newton.raphson.bic(data=df[2:(last.age+1),2017-1961+1],offset=log(Ef[2:(last.age+1),2017-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(1),lambda=c(2,8),typsize=c(1,1),fscale=1,stepmax=5)

m.2010.100<-newton.raphson.bic(data=dm[2:(last.age+1),2010-1961+1],offset=log(Em[2:(last.age+1),2010-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(0),lambda=c(4),typsize=c(1),fscale=1,stepmax=10)
m.2011.100<-newton.raphson.bic(data=dm[2:(last.age+1),2011-1961+1],offset=log(Em[2:(last.age+1),2011-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(0),lambda=c(4),typsize=c(1),fscale=1,stepmax=10)
m.2012.100<-newton.raphson.bic(data=dm[2:(last.age+1),2012-1961+1],offset=log(Em[2:(last.age+1),2012-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(0),lambda=c(4),typsize=c(1),fscale=1,stepmax=10)
m.2010.expo.100<-newton.raphson.bic(data=dm[2:(last.age+1),2010-1961+1],offset=log(Em[2:(last.age+1),2010-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(1),lambda=c(2,8),typsize=c(1,1),fscale=1,stepmax=5)
m.2011.expo.100<-newton.raphson.bic(data=dm[2:(last.age+1),2011-1961+1],offset=log(Em[2:(last.age+1),2011-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(1),lambda=c(2,8),typsize=c(1,1),fscale=1,stepmax=5)
m.2012.expo.100<-newton.raphson.bic(data=dm[2:(last.age+1),2012-1961+1],offset=log(Em[2:(last.age+1),2012-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(1),lambda=c(2,8),typsize=c(1,1),fscale=1,stepmax=5)
f.2010.100<-newton.raphson.bic(data=df[2:(last.age+1),2010-1961+1],offset=log(Ef[2:(last.age+1),2010-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(0),lambda=c(2),typsize=c(1),fscale=1,stepmax=2)
f.2011.100<-newton.raphson.bic(data=df[2:(last.age+1),2011-1961+1],offset=log(Ef[2:(last.age+1),2011-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(0),lambda=c(4),typsize=c(1),fscale=1,stepmax=10)
f.2012.100<-newton.raphson.bic(data=df[2:(last.age+1),2012-1961+1],offset=log(Ef[2:(last.age+1),2012-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(0),lambda=c(4),typsize=c(1),fscale=1,stepmax=10)
f.2010.expo.100<-newton.raphson.bic(data=df[2:(last.age+1),2010-1961+1],offset=log(Ef[2:(last.age+1),2010-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(1),lambda=c(2,8),typsize=c(1,1),fscale=1,stepmax=5)
f.2011.expo.100<-newton.raphson.bic(data=df[2:(last.age+1),2011-1961+1],offset=log(Ef[2:(last.age+1),2011-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(1),lambda=c(2,8),typsize=c(1,1),fscale=1,stepmax=5)
f.2012.expo.100<-newton.raphson.bic(data=df[2:(last.age+1),2012-1961+1],offset=log(Ef[2:(last.age+1),2012-1961+1]),DXX=A,P.list=list(list(P)),P.index=list(1),exp.index=c(1),lambda=c(2,8),typsize=c(1,1),fscale=1,stepmax=5)

j.2007<-newton.raphson.bic(data=c(dm[2:(last.age+1),2007-1961+1],df[2:(last.age+1),2007-1961+1]),offset=log(c(Em[2:(last.age+1),2007-1961+1],Ef[2:(last.age+1),2007-1961+1])),DXX=adiag(A,A),P.list=list(list(P),list(P),list(PD[-(1:8),])),P.index=list(1,no.basis+1,1),exp.index=c(1,1,1),lambda=c(2,6,2,6,2,6),typsize=c(1,1,1,1,1,1),fscale=1,stepmax=50)
j.2010<-newton.raphson.bic(data=c(dm[2:(last.age+1),2010-1961+1],df[2:(last.age+1),2010-1961+1]),offset=log(c(Em[2:(last.age+1),2010-1961+1],Ef[2:(last.age+1),2010-1961+1])),DXX=adiag(A,A),P.list=list(list(P),list(P),list(PD[-(1:8),])),P.index=list(1,no.basis+1,1),exp.index=c(1,1,1),lambda=c(2,6,2,6,2,6),typsize=c(1,1,1,1,1,1),fscale=1,stepmax=50)
j.2011<-newton.raphson.bic(data=c(dm[2:(last.age+1),2011-1961+1],df[2:(last.age+1),2011-1961+1]),offset=log(c(Em[2:(last.age+1),2011-1961+1],Ef[2:(last.age+1),2011-1961+1])),DXX=adiag(A,A),P.list=list(list(P),list(P),list(PD[-(1:8),])),P.index=list(1,no.basis+1,1),exp.index=c(1,1,1),lambda=c(2,6,2,6,2,6),typsize=c(1,1,1,1,1,1),fscale=1,stepmax=50)
j.2012<-newton.raphson.bic(data=c(dm[2:(last.age+1),2012-1961+1],df[2:(last.age+1),2012-1961+1]),offset=log(c(Em[2:(last.age+1),2012-1961+1],Ef[2:(last.age+1),2012-1961+1])),DXX=adiag(A,A),P.list=list(list(P),list(P),list(PD[-(1:8),])),P.index=list(1,no.basis+1,1),exp.index=c(1,1,1),lambda=c(2,6,2,6,2,6),typsize=c(1,1,1,1,1,1),fscale=1,stepmax=50)
j.2017<-newton.raphson.bic(data=c(dm[2:(last.age+1),2017-1961+1],df[2:(last.age+1),2017-1961+1]),offset=log(c(Em[2:(last.age+1),2017-1961+1],Ef[2:(last.age+1),2017-1961+1])),DXX=adiag(A,A),P.list=list(list(P),list(P),list(PD[-(1:8),])),P.index=list(1,no.basis+1,1),exp.index=c(1,1,1),lambda=c(2,6,2,6,2,6),typsize=c(1,1,1,1,1,1),fscale=1,stepmax=100)

graduation.fit<-function(x,exprate.index,data,offset){
lambda<-x$sp
if(exprate.index==1) lambda.weights<-exp(lambda[1]+lambda[2]*seq(0,1,length=nrow(P))) else lambda.weights<-exp(lambda[1])
B<-as.numeric(lambda.weights^0.5)*P
fit<-IRLS(data=data,DXX=A,B=B,offset=offset)
fit
}

graduation.fit.joint<-function(x,data,offset){
lambda<-x$sp
lambda.weightsm<-exp(lambda[1]+lambda[2]*seq(0,1,length=nrow(P)))
lambda.weightsf<-exp(lambda[3]+lambda[4]*seq(0,1,length=nrow(P)))
lambda.weightsd<-exp(lambda[5]+lambda[6]*seq(0,1,length=nrow(PD[-(1:8),])))
#lambda.weightsm<-c(exp(lambda[1]+lambda[2]*seq(0,1,length=no.basis-2)),rep(exp(lambda[1]+lambda[2]),extra.knots))
#lambda.weightsf<-c(exp(lambda[3]+lambda[4]*seq(0,1,length=no.basis-2)),rep(exp(lambda[3]+lambda[4]),extra.knots))
#lambda.weightsd<-c(exp(lambda[5]+lambda[6]*seq(0,1,length=no.basis-8)),rep(exp(lambda[5]+lambda[6]),extra.knots))
Bm<-as.numeric(lambda.weightsm^0.5)*P
Bf<-as.numeric(lambda.weightsf^0.5)*P
Bd<-as.numeric(lambda.weightsd^0.5)*PD[-(1:8),]
B<-rbind(adiag(Bm,Bf),Bd)
fit<-IRLS(y=data,X=adiag(A,A),B=B,offset=offset)
fit
}

extra.knots=6
knots<-seq(0,1,length=no.basis-2)
dk<-knots[2]-knots[1]	
knots<-c(knots[1]-dk*(3:1),knots,knots[no.basis-2]+dk*(1:(3+extra.knots)))
age<-seq(0,1,length=last.age)

A<-c()
for(j in 1:(no.basis+extra.knots)) {
A<-cbind(A,bspline(age,knots,j))
}
P<-diff(diag(no.basis+extra.knots),differences=2)
PD<-t(c(1,-1))%x%diag(no.basis+extra.knots)

m.1980.extra<-graduation.fit(m.1980,0,data=dm[2:(last.age+1),1980-1961+1],offset=log(Em[2:(last.age+1),1980-1961+1]))
m.2010.extra<-graduation.fit(m.2010,0,data=dm[2:(last.age+1),2010-1961+1],offset=log(Em[2:(last.age+1),2010-1961+1]))
m.2011.extra<-graduation.fit(m.2011,0,data=dm[2:(last.age+1),2011-1961+1],offset=log(Em[2:(last.age+1),2011-1961+1]))
m.2012.extra<-graduation.fit(m.2012,0,data=dm[2:(last.age+1),2012-1961+1],offset=log(Em[2:(last.age+1),2012-1961+1]))
m.2017.extra<-graduation.fit(m.2017,0,data=dm[2:(last.age+1),2017-1961+1],offset=log(Em[2:(last.age+1),2017-1961+1]))
m.1980.extra.expo<-graduation.fit(m.1980.expo,1,data=dm[2:(last.age+1),1980-1961+1],offset=log(Em[2:(last.age+1),1980-1961+1]))
m.2007.extra.expo<-graduation.fit(m.2007.expo,1,data=dm[2:(last.age+1),2007-1961+1],offset=log(Em[2:(last.age+1),2007-1961+1]))
m.2010.extra.expo<-graduation.fit(m.2010.expo,1,data=dm[2:(last.age+1),2010-1961+1],offset=log(Em[2:(last.age+1),2010-1961+1]))
m.2011.extra.expo<-graduation.fit(m.2011.expo,1,data=dm[2:(last.age+1),2011-1961+1],offset=log(Em[2:(last.age+1),2011-1961+1]))
m.2012.extra.expo<-graduation.fit(m.2012.expo,1,data=dm[2:(last.age+1),2012-1961+1],offset=log(Em[2:(last.age+1),2012-1961+1]))
m.2017.extra.expo<-graduation.fit(m.2017.expo,1,data=dm[2:(last.age+1),2017-1961+1],offset=log(Em[2:(last.age+1),2017-1961+1]))
f.2010.extra<-graduation.fit(f.2010,0,data=df[2:(last.age+1),2010-1961+1],offset=log(Ef[2:(last.age+1),2010-1961+1]))
f.2011.extra<-graduation.fit(f.2011,0,data=df[2:(last.age+1),2011-1961+1],offset=log(Ef[2:(last.age+1),2011-1961+1]))
f.2012.extra<-graduation.fit(f.2012,0,data=df[2:(last.age+1),2012-1961+1],offset=log(Ef[2:(last.age+1),2012-1961+1]))
f.1980.extra.expo<-graduation.fit(f.1980.expo,1,data=df[2:(last.age+1),1980-1961+1],offset=log(Ef[2:(last.age+1),1980-1961+1]))
f.2007.extra.expo<-graduation.fit(f.2007.expo,1,data=df[2:(last.age+1),2007-1961+1],offset=log(Ef[2:(last.age+1),2007-1961+1]))
f.2010.extra.expo<-graduation.fit(f.2010.expo,1,data=df[2:(last.age+1),2010-1961+1],offset=log(Ef[2:(last.age+1),2010-1961+1]))
f.2011.extra.expo<-graduation.fit(f.2011.expo,1,data=df[2:(last.age+1),2011-1961+1],offset=log(Ef[2:(last.age+1),2011-1961+1]))
f.2012.extra.expo<-graduation.fit(f.2012.expo,1,data=df[2:(last.age+1),2012-1961+1],offset=log(Ef[2:(last.age+1),2012-1961+1]))
f.2017.extra.expo<-graduation.fit(f.2017.expo,1,data=df[2:(last.age+1),2017-1961+1],offset=log(Ef[2:(last.age+1),2017-1961+1]))

m.2010.extra.100<-graduation.fit(m.2010,0,data=dm[2:(last.age+1),2010-1961+1],offset=log(Em[2:(last.age+1),2010-1961+1]))
m.2011.extra.100<-graduation.fit(m.2011,0,data=dm[2:(last.age+1),2011-1961+1],offset=log(Em[2:(last.age+1),2011-1961+1]))
m.2012.extra.100<-graduation.fit(m.2012,0,data=dm[2:(last.age+1),2012-1961+1],offset=log(Em[2:(last.age+1),2012-1961+1]))
m.2010.extra.expo.100<-graduation.fit(m.2010.expo,1,data=dm[2:(last.age+1),2010-1961+1],offset=log(Em[2:(last.age+1),2010-1961+1]))
m.2011.extra.expo.100<-graduation.fit(m.2011.expo,1,data=dm[2:(last.age+1),2011-1961+1],offset=log(Em[2:(last.age+1),2011-1961+1]))
m.2012.extra.expo.100<-graduation.fit(m.2012.expo,1,data=dm[2:(last.age+1),2012-1961+1],offset=log(Em[2:(last.age+1),2012-1961+1]))
f.2010.extra.100<-graduation.fit(f.2010,0,data=df[2:(last.age+1),2010-1961+1],offset=log(Ef[2:(last.age+1),2010-1961+1]))
f.2011.extra.100<-graduation.fit(f.2011,0,data=df[2:(last.age+1),2011-1961+1],offset=log(Ef[2:(last.age+1),2011-1961+1]))
f.2012.extra.100<-graduation.fit(f.2012,0,data=df[2:(last.age+1),2012-1961+1],offset=log(Ef[2:(last.age+1),2012-1961+1]))
f.2010.extra.expo.100<-graduation.fit(f.2010.expo,1,data=df[2:(last.age+1),2010-1961+1],offset=log(Ef[2:(last.age+1),2010-1961+1]))
f.2011.extra.expo.100<-graduation.fit(f.2011.expo,1,data=df[2:(last.age+1),2011-1961+1],offset=log(Ef[2:(last.age+1),2011-1961+1]))
f.2012.extra.expo.100<-graduation.fit(f.2012.expo,1,data=df[2:(last.age+1),2012-1961+1],offset=log(Ef[2:(last.age+1),2012-1961+1]))

j.2007.extra<-graduation.fit.joint(j.2007,data=c(dm[2:(last.age+1),2007-1961+1],df[2:(last.age+1),2007-1961+1]),offset=log(c(Em[2:(last.age+1),2007-1961+1],Ef[2:(last.age+1),2007-1961+1])))
j.2010.extra<-graduation.fit.joint(j.2010,data=c(dm[2:(last.age+1),2010-1961+1],df[2:(last.age+1),2010-1961+1]),offset=log(c(Em[2:(last.age+1),2010-1961+1],Ef[2:(last.age+1),2010-1961+1])))
j.2011.extra<-graduation.fit.joint(j.2011,data=c(dm[2:(last.age+1),2011-1961+1],df[2:(last.age+1),2011-1961+1]),offset=log(c(Em[2:(last.age+1),2011-1961+1],Ef[2:(last.age+1),2011-1961+1])))
j.2012.extra<-graduation.fit.joint(j.2012,data=c(dm[2:(last.age+1),2012-1961+1],df[2:(last.age+1),2012-1961+1]),offset=log(c(Em[2:(last.age+1),2012-1961+1],Ef[2:(last.age+1),2012-1961+1])))
j.2017.extra<-graduation.fit.joint(j.2017,data=c(dm[2:(last.age+1),2017-1961+1],df[2:(last.age+1),2017-1961+1]),offset=log(c(Em[2:(last.age+1),2017-1961+1],Ef[2:(last.age+1),2017-1961+1])))

A.extra<-c()
for(j in 1:(no.basis+extra.knots)) {
A.extra<-cbind(A.extra,bspline(0:119/(last.age-1),knots,j))
}

A.extra.100<-c()
for(j in 1:(no.basis+extra.knots)) {
A.extra.100<-cbind(A.extra.100,bspline(0:119/(100-1),knots,j))
}

par(mfrow=c(1,3),cex.lab=2,cex.axis=1.5,cex.main=1.5,mar=c(5,5,4,2))
plot(A.extra.100%*%m.2010.extra.100$b,type="l",col=4,lty=2,main="Males 2010",ylab="Log Mortality Rates",xlab="Age",ylim=c(-9,1))
lines(A.extra%*%m.2010.extra$b,type="l",col=4)
points(log(dm/Em)[2:106,2010-1961+1],col=4)
points(101:105,log(dm/Em)[102:106,2010-1961+1],col=2)
legend("topleft",c("Age 1 to 100","Age 1 to 105"),lty=c(2,1),col=c(4,4),bty="n",x.intersp=0.3,inset=c(-0.08,0))
plot(A.extra.100%*%m.2011.extra.100$b,type="l",col=4,lty=2,main="Males 2011",ylab="Log Mortality Rates",xlab="Age",ylim=c(-9,1))
lines(A.extra%*%m.2011.extra$b,type="l",col=4)
points(log(dm/Em)[2:106,2011-1961+1],col=4)
points(101:105,log(dm/Em)[102:106,2011-1961+1],col=2)
legend("topleft",c("Age 1 to 100","Age 1 to 105"),lty=c(2,1),col=c(4,4),bty="n",x.intersp=0.3,inset=c(-0.08,0))
plot(A.extra.100%*%m.2012.extra.100$b,type="l",col=4,lty=2,main="Males 2012",ylab="Log Mortality Rates",xlab="Age",ylim=c(-9,1))
lines(A.extra%*%m.2012.extra$b,type="l",col=4)
points(log(dm/Em)[2:106,2012-1961+1],col=4)
points(101:105,log(dm/Em)[102:106,2012-1961+1],col=2)
legend("topleft",c("Age 1 to 100","Age 1 to 105"),lty=c(2,1),col=c(4,4),bty="n",x.intersp=0.3,inset=c(-0.08,0))

plot(A.extra.100%*%f.2010.extra.100$b,type="l",col=2,lty=2,main="Females 2010",ylab="Log Mortality Rates",xlab="Age",ylim=c(-9.5,1))
lines(A.extra%*%f.2010.extra$b,type="l",col=2)
points(log(df/Ef)[2:106,2010-1961+1],col=2)
points(101:105,log(df/Ef)[102:106,2010-1961+1],col=4)
legend("topleft",c("Age 1 to 100","Age 1 to 105"),lty=c(2,1),col=c(2,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))
plot(A.extra.100%*%f.2011.extra.100$b,type="l",col=2,lty=2,main="Females 2011",ylab="Log Mortality Rates",xlab="Age",ylim=c(-9.5,1))
lines(A.extra%*%f.2011.extra$b,type="l",col=2)
points(log(df/Ef)[2:106,2011-1961+1],col=2)
points(101:105,log(df/Ef)[102:106,2011-1961+1],col=4)
legend("topleft",c("Age 1 to 100","Age 1 to 105"),lty=c(2,1),col=c(2,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))
plot(A.extra.100%*%f.2012.extra.100$b,type="l",col=2,lty=2,main="Females 2012",ylab="Log Mortality Rates",xlab="Age",ylim=c(-9.5,1))
lines(A.extra%*%f.2012.extra$b,type="l",col=2)
points(log(df/Ef)[2:106,2012-1961+1],col=2)
points(101:105,log(df/Ef)[102:106,2012-1961+1],col=4)
legend("topleft",c("Age 1 to 100","Age 1 to 105"),lty=c(2,1),col=c(2,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))

par(mfrow=c(1,3),cex.lab=2,cex.axis=1.5,cex.main=1.5,mar=c(5,5,4,2))
plot(A.extra%*%m.2010.extra$b,type="l",col=4,ylim=c(-9.5,1),main="Year 2010",ylab="Log Mortality Rates",xlab="Age")
lines(A.extra%*%f.2010.extra$b,type="l",col=2)
legend("topleft",c("Males","Females"),lty=c(1,1),col=c(4,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))
plot(A.extra%*%m.2011.extra$b,type="l",col=4,ylim=c(-9.5,1),main="Year 2011",ylab="Log Mortality Rates",xlab="Age")
lines(A.extra%*%f.2011.extra$b,type="l",col=2)
legend("topleft",c("Males","Females"),lty=c(1,1),col=c(4,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))
plot(A.extra%*%m.2012.extra$b,type="l",col=4,ylim=c(-9.5,1),main="Year 2012",ylab="Log Mortality Rates",xlab="Age")
lines(A.extra%*%f.2012.extra$b,type="l",col=2)
legend("topleft",c("Males","Females"),lty=c(1,1),col=c(4,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))

par(mfrow=c(1,3),cex.lab=2,cex.axis=1.5,cex.main=1.5,mar=c(5,5,4,2))
plot(A.extra%*%m.2010.extra.expo$b,type="l",col=4,ylim=c(-9.5,1),main="Year 2010",ylab="Log Mortality Rates",xlab="Age")
lines(A.extra%*%f.2010.extra.expo$b,type="l",col=2)
lines(A.extra%*%m.2010.extra$b,type="l",col=4,lty=2)
lines(A.extra%*%f.2010.extra$b,type="l",col=2,lty=2)
legend("topleft",c("Males (local penalty)","Females (local penalty)","Males (global penalty)","Females (global penalty)"),lty=c(1,1,2,2),col=c(4,2,4,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))
plot(A.extra%*%m.2011.extra.expo$b,type="l",col=4,ylim=c(-9.5,1),main="Year 2011",ylab="Log Mortality Rates",xlab="Age")
lines(A.extra%*%f.2011.extra.expo$b,type="l",col=2)
lines(A.extra%*%m.2011.extra$b,type="l",col=4,lty=2)
lines(A.extra%*%f.2011.extra$b,type="l",col=2,lty=2)
legend("topleft",c("Males (local penalty)","Females (local penalty)","Males (global penalty)","Females (global penalty)"),lty=c(1,1,2,2),col=c(4,2,4,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))
plot(A.extra%*%m.2012.extra.expo$b,type="l",col=4,ylim=c(-9.5,1),main="Year 2012",ylab="Log Mortality Rates",xlab="Age")
lines(A.extra%*%f.2012.extra.expo$b,type="l",col=2)
lines(A.extra%*%m.2012.extra$b,type="l",col=4,lty=2)
lines(A.extra%*%f.2012.extra$b,type="l",col=2,lty=2)
legend("topleft",c("Males (local penalty)","Females (local penalty)","Males (global penalty)","Females (global penalty)"),lty=c(1,1,2,2),col=c(4,2,4,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))

plot(A.extra%*%m.1980.extra$b,type="l",col=4,ylim=c(-8.5,0.5),main="Males 1980",ylab="Log Mortality Rates",xlab="Age",lty=2)
lines(A.extra%*%m.1980.extra.expo$b,type="l",col=4)
legend("topleft",c("Males (local penalty)","Males (global penalty)"),lty=c(1,2),col=c(4,4),bty="n",x.intersp=0.3,inset=c(-0.008,0))

par(mfrow=c(1,2),cex.lab=2,cex.axis=1.5,cex.main=1.5,mar=c(5,5,4,2))
plot(A.extra%*%m.2007.extra.expo$b,type="l",col=4,ylim=c(-9,1),main="Year 2007",ylab="Log Mortality Rates",xlab="Age")
lines(A.extra%*%f.2007.extra.expo$b,type="l",col=2)
legend("topleft",c("Males","Females"),lty=c(1,1),col=c(4,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))
plot(A.extra%*%m.2017.extra.expo$b,type="l",col=4,ylim=c(-9.5,0.5),main="Year 2017",ylab="Log Mortality Rates",xlab="Age")
lines(A.extra%*%f.2017.extra.expo$b,type="l",col=2)
legend("topleft",c("Males","Females"),lty=c(1,1),col=c(4,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))

par(mfrow=c(1,3),cex.lab=2,cex.axis=1.5,cex.main=1.5,mar=c(5,5,4,2))
plot(A.extra%*%j.2010.extra$b[1:(no.basis+extra.knots)],type="l",col=4,ylim=c(-9.5,1),main="Year 2010",ylab="Log Mortality Rates",xlab="Age")
lines(A.extra%*%j.2010.extra$b[no.basis+extra.knots+1:(no.basis+extra.knots)],type="l",col=2)
lines(A.extra%*%m.2010.extra.expo$b,type="l",col=4,lty=2)
lines(A.extra%*%f.2010.extra.expo$b,type="l",col=2,lty=2)
legend("topleft",c("Males (joint)","Females (joint)","Males (separate)","Females (separate)"),lty=c(1,1,2,2),col=c(4,2,4,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))
plot(A.extra%*%j.2011.extra$b[1:(no.basis+extra.knots)],type="l",col=4,ylim=c(-9.5,1),main="Year 2011",ylab="Log Mortality Rates",xlab="Age")
lines(A.extra%*%j.2011.extra$b[no.basis+extra.knots+1:(no.basis+extra.knots)],type="l",col=2)
lines(A.extra%*%m.2011.extra.expo$b,type="l",col=4,lty=2)
lines(A.extra%*%f.2011.extra.expo$b,type="l",col=2,lty=2)
legend("topleft",c("Males (joint)","Females (joint)","Males (separate)","Females (separate)"),lty=c(1,1,2,2),col=c(4,2,4,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))
plot(A.extra%*%j.2012.extra$b[1:(no.basis+extra.knots)],type="l",col=4,ylim=c(-9.5,1),main="Year 2012",ylab="Log Mortality Rates",xlab="Age")
lines(A.extra%*%j.2012.extra$b[no.basis+extra.knots+1:(no.basis+extra.knots)],type="l",col=2)
lines(A.extra%*%m.2012.extra.expo$b,type="l",col=4,lty=2)
lines(A.extra%*%f.2012.extra.expo$b,type="l",col=2,lty=2)
legend("topleft",c("Males (joint)","Females (joint)","Males (separate)","Females (separate)"),lty=c(1,1,2,2),col=c(4,2,4,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))

par(mfrow=c(1,2),cex.lab=2,cex.axis=1.5,cex.main=1.5,mar=c(5,5,4,2))
plot(A.extra%*%j.2007.extra$b[1:(no.basis+extra.knots)],type="l",col=4,ylim=c(-9.5,1),main="Year 2007",ylab="Log Mortality Rates",xlab="Age")
lines(A.extra%*%j.2007.extra$b[no.basis+extra.knots+1:(no.basis+extra.knots)],type="l",col=2)
lines(A.extra%*%m.2007.extra.expo$b,type="l",col=4,lty=2)
lines(A.extra%*%f.2007.extra.expo$b,type="l",col=2,lty=2)
legend("topleft",c("Males (joint)","Females (joint)","Males (separate)","Females (separate)"),lty=c(1,1,2,2),col=c(4,2,4,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))
plot(A.extra%*%j.2017.extra$b[1:(no.basis+extra.knots)],type="l",col=4,ylim=c(-9.5,1),main="Year 2017",ylab="Log Mortality Rates",xlab="Age")
lines(A.extra%*%j.2017.extra$b[no.basis+extra.knots+1:(no.basis+extra.knots)],type="l",col=2)
lines(A.extra%*%m.2017.extra.expo$b,type="l",col=4,lty=2)
lines(A.extra%*%f.2017.extra.expo$b,type="l",col=2,lty=2)
legend("topleft",c("Males (joint)","Females (joint)","Males (separate)","Females (separate)"),lty=c(1,1,2,2),col=c(4,2,4,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))


par(mfrow=c(1,2),cex.lab=2,cex.axis=1.5,cex.main=1.5,mar=c(5,5,4,2))
plot(A.extra%*%j.2007.extra$b[1:(no.basis+extra.knots)],type="l",col=4,ylim=c(-9.5,1),main="Year 2007",ylab="Log Mortality Rates",xlab="Age")
lines(A.extra%*%j.2007.extra$b[no.basis+extra.knots+1:(no.basis+extra.knots)],type="l",col=2)
legend("topleft",c("Males (joint)","Females (joint)"),lty=c(1,1),col=c(4,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))
plot(c(1:7,9,9:46),c(rep(0,8),exp(j.2007$sp[5]+j.2007$sp[6]*seq(0,1,length=nrow(PD[-(1:8),])))),type="l",xlab="Basis",ylab="Penalty",main=expression(zeta^D~(i)))

#######NNLS
library(lsei)
P.m<-cbind(diff(diag(no.basis+extra.knots),differences=2),matrix(0,no.basis+extra.knots-2,no.basis+extra.knots))
P.f<-cbind(matrix(0,no.basis+extra.knots-2,no.basis+extra.knots),diff(diag(no.basis+extra.knots),differences=2))
P.d<-t(c(1,-1))%x%diag(no.basis+extra.knots)[-(1:8),]

para.trans<-rbind(cbind(matrix(0,no.basis+extra.knots,no.basis+extra.knots),diag(no.basis+extra.knots)),t(c(1,-1))%x%diag(no.basis+extra.knots))
inv.para.trans<-solve(para.trans)
DX.perm<-adiag(A,A)%*%inv.para.trans
P.m.perm<-P.m%*%inv.para.trans
P.f.perm<-P.f%*%inv.para.trans
P.d.perm<-P.d%*%inv.para.trans

IRLS.nnls<-function (y,X,B,mu,offset) {
mu<-ifelse(mu==0,mu+1e-3,mu)
dev<-ll.sat<-sum(dpois(y,y,log=TRUE))
eta<-log(mu)
converged<-FALSE
while(!converged) {
	z<-(y-mu)/mu+eta-offset
	w<-mu
	fit<-pnnls(rbind(as.numeric(w^0.5)*X,B),c(as.numeric(w^0.5)*z,rep(0,nrow(B))),k=no.basis+extra.knots)
	eta<-X%*%fit$x+offset
	mu<-exp(eta)
	old.dev<-dev
	dev<-2*(ll.sat-sum(dpois(y,mu,log=TRUE)))
	if(abs(dev-old.dev)<1e-6*dev) converged<-TRUE
}
b=fit$x
fv=mu
list(b=b,fv=fv,dev=dev)
}

extrapolate.nnls<-function(x,data,offset){
lambda<-x$sp
weightsm<-exp(lambda[1]+lambda[2]*seq(0,1,length=no.basis+extra.knots-2))
weightsf<-exp(lambda[3]+lambda[4]*seq(0,1,length=no.basis+extra.knots-2))
weightsd<-exp(lambda[5]+lambda[6]*seq(0,1,length=nrow(P.d)))
B1<-sqrt(weightsm)*P.m.perm
B2<-sqrt(weightsf)*P.f.perm
B3<-sqrt(weightsd)*P.d.perm
B<-rbind(B1,B2,B3)
fit<-IRLS.nnls(data,X=DX.perm,B=B,mu=data,offset=offset)
fit
}

j.2010.nnls.extra<-extrapolate.nnls(j.2010,data=c(dm[2:(last.age+1),2010-1961+1],df[2:(last.age+1),2010-1961+1]),offset=log(c(Em[2:(last.age+1),2010-1961+1],Ef[2:(last.age+1),2010-1961+1])))
j.2011.nnls.extra<-extrapolate.nnls(j.2011,data=c(dm[2:(last.age+1),2011-1961+1],df[2:(last.age+1),2011-1961+1]),offset=log(c(Em[2:(last.age+1),2011-1961+1],Ef[2:(last.age+1),2011-1961+1])))
j.2012.nnls.extra<-extrapolate.nnls(j.2012,data=c(dm[2:(last.age+1),2012-1961+1],df[2:(last.age+1),2012-1961+1]),offset=log(c(Em[2:(last.age+1),2012-1961+1],Ef[2:(last.age+1),2012-1961+1])))

par(mfrow=c(1,3),cex.lab=2,cex.axis=1.5,cex.main=1.5,mar=c(5,5,4,2))
plot(A.extra%*%(inv.para.trans%*%j.2010.nnls.extra$b)[1:(no.basis+extra.knots)],type="l",col=4,ylim=c(-9.5,1),main="Year 2010",ylab="Log Mortality Rates",xlab="Age")
lines(A.extra%*%(inv.para.trans%*%j.2010.nnls.extra$b)[no.basis+extra.knots+1:(no.basis+extra.knots)],type="l",col=2)
legend("topleft",c("Males","Females"),lty=c(1,1),col=c(4,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))
plot(A.extra%*%(inv.para.trans%*%j.2011.nnls.extra$b)[1:(no.basis+extra.knots)],type="l",col=4,ylim=c(-9.5,1),main="Year 2011",ylab="Log Mortality Rates",xlab="Age")
lines(A.extra%*%(inv.para.trans%*%j.2011.nnls.extra$b)[no.basis+extra.knots+1:(no.basis+extra.knots)],type="l",col=2)
legend("topleft",c("Males","Females"),lty=c(1,1),col=c(4,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))
plot(A.extra%*%(inv.para.trans%*%j.2012.nnls.extra$b)[1:(no.basis+extra.knots)],type="l",col=4,ylim=c(-9.5,1),main="Year 2012",ylab="Log Mortality Rates",xlab="Age")
lines(A.extra%*%(inv.para.trans%*%j.2012.nnls.extra$b)[no.basis+extra.knots+1:(no.basis+extra.knots)],type="l",col=2)
legend("topleft",c("Males","Females"),lty=c(1,1),col=c(4,2),bty="n",x.intersp=0.3,inset=c(-0.08,0))

par(mfrow=c(1,2))
plot(A.extra%*%m.2011.extra$b,ylab="Log Mortality Rates",xlab="Age",main="England and Wales Males 2011 \n (on log scale)",col=4,cex.lab=1.5,cex.main=1.5,type="l")
points(log(dm/Em)[2:106,2011-1961+1],col=4)
plot(A.extra%*%f.2011.extra$b,ylab="Log Mortality Rates",xlab="Age",main="England and Wales Females 2011 \n (on log scale)",col=2,cex.lab=1.5,cex.main=1.5,type="l")
points(log(df/Ef)[2:106,2011-1961+1],col=2)

plot(A.extra%*%m.2011.extra.expo$b,ylab="Log Mortality Rates",xlab="Age",main="England and Wales Males 2011 \n (on log scale)",col=4,cex.lab=1.5,cex.main=1.5,type="l")
points(log(dm/Em)[2:106,2011-1961+1],col=4)
plot(exp(m.2011.expo$sp[1]+m.2011.expo$sp[2]*seq(0,1,length=38)),type="l",xlab="Basis",ylab="Penalty",main=expression(zeta (i)),cex.lab=1.5,cex.main=1.5)
