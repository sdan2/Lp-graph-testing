library(igraph)
library(Matrix)
library(tikzDevice) 


sim_Z_ER <- function(n = 100, p=0.5, m = 32, r = 4)
{
  c = m/(2*r)
  b = n/2
  Prod_var = as(matrix(1, n, n), "sparseMatrix")
  Prod_stat = as(matrix(1, n, n), "sparseMatrix")
  Diff = as(matrix(0, n, n), "sparseMatrix")
  for(l in 1:(2*r))
  {
    Temp = Diff
    Diff = as(matrix(0, n, n), "sparseMatrix")
    Sum = as(matrix(0, n, n), "sparseMatrix")
    for(s in 1:c)
    {
      A = as_adj(sample_gnp(n, p))
      B = as_adj(sample_gnp(n, p))
      Diff = Diff + as(Matrix(A - B, sparse = F), "sparseMatrix")
      Sum = Sum + Matrix(A + B, sparse = T)
    }
    if(l%%2 == 0)
    {
      Prod_stat = Prod_stat*(as(Matrix(Temp + Diff, sparse = F), "sparseMatrix"))
    }
    if(l <= r)
    {
      Prod_var = Prod_var*(Sum/(2*c))
    }
    else
    {
      Prod_var = Prod_var*(1 - Sum/(2*c))
    }
  }
  T_stat = 0
  S_stat = 0
  for(i in 1:(n - 1))
  {
    for(j in (i + 1):n)
    {
      T_stat = T_stat + Prod_stat[i,j]
      S_stat = S_stat + Prod_var[i,j]
    }
  }
  S_stat = S_stat*(4*c)^r
  Z_stat = T_stat/sqrt(S_stat)
  return(Z_stat)
}

sim_Z_PB <- function(n = 100, a=8, b=2, m = 32, r = 4)
{
  c = m/(2*r)
  b = n/2
  Prod_var = as(matrix(1, n, n), "sparseMatrix")
  Prod_stat = as(matrix(1, n, n), "sparseMatrix")
  Diff = as(matrix(0, n, n), "sparseMatrix")
  for(l in 1:(2*r))
  {
    Temp = Diff
    Diff = as(matrix(0, n, n), "sparseMatrix")
    Sum = as(matrix(0, n, n), "sparseMatrix")
    for(s in 1:c)
    {
      A = as_adj(sample_sbm(n, pref.matrix = matrix(c(a/n, b/n, b/n, a/n), 2, 2), block.sizes = c(n/2, n/2)))
      B = as_adj(sample_sbm(n, pref.matrix = matrix(c(a/n, b/n, b/n, a/n), 2, 2), block.sizes = c(n/2, n/2)))
      Diff = Diff + as(Matrix(A - B, sparse = F), "sparseMatrix")
      Sum = Sum + Matrix(A + B, sparse = T)
    }
    if(l%%2 == 0)
    {
      Prod_stat = Prod_stat*(as(Matrix(Temp + Diff, sparse = F), "sparseMatrix"))
    }
    if(l <= r)
    {
      Prod_var = Prod_var*(Sum/(2*c))
    }
    else
    {
      Prod_var = Prod_var*(1 - Sum/(2*c))
    }
  }
  T_stat = 0
  S_stat = 0
  for(i in 1:(n - 1))
  {
    for(j in (i + 1):n)
    {
      T_stat = T_stat + Prod_stat[i,j]
      S_stat = S_stat + Prod_var[i,j]
    }
  }
  S_stat = S_stat*(4*c)^r
  Z_stat = T_stat/sqrt(S_stat)
  return(Z_stat)
}

graphon <- function(x, y){exp(x + y)/(1 + exp(x + y))}

sample_beta <- function(n = 100, rad = 1)
{
  beta = rnorm(n)
  beta = rad*beta/norm(beta, type = '2')
  adj = matrix(0, n, n)
  for(.i in 2:n)
  {                         
    for(.j in 1:(.i-1))
    {
      .p = graphon(beta[.i], beta[.j])
      adj[.i, .j] = rbinom(1, 1, .p)
    }
  }
  adjsymm = as(adj + t(adj), 'sparseMatrix')
  return(adjsymm)
}

sim_Z_beta <- function(n = 100, rad = 1, m = 32, r = 4)
{
  c = m/(2*r)
  b = n/2
  Prod_var = as(matrix(1, n, n), "sparseMatrix")
  Prod_stat = as(matrix(1, n, n), "sparseMatrix")
  Diff = as(matrix(0, n, n), "sparseMatrix")
  for(l in 1:(2*r))
  {
    Temp = Diff
    Diff = as(matrix(0, n, n), "sparseMatrix")
    Sum = as(matrix(0, n, n), "sparseMatrix")
    for(s in 1:c)
    {
      A = sample_beta(n, rad)
      B = sample_beta(n, rad)
      Diff = Diff + as(Matrix(A - B, sparse = F), "sparseMatrix")
      Sum = Sum + Matrix(A + B, sparse = T)
    }
    if(l%%2 == 0)
    {
      Prod_stat = Prod_stat*(as(Matrix(Temp + Diff, sparse = F), "sparseMatrix"))
    }
    if(l <= r)
    {
      Prod_var = Prod_var*(Sum/(2*c))
    }
    else
    {
      Prod_var = Prod_var*(1 - Sum/(2*c))
    }
  }
  T_stat = 0
  S_stat = 0
  for(i in 1:(n - 1))
  {
    for(j in (i + 1):n)
    {
      T_stat = T_stat + Prod_stat[i,j]
      S_stat = S_stat + Prod_var[i,j]
    }
  }
  S_stat = S_stat*(4*c)^r
  Z_stat = T_stat/sqrt(S_stat)
  return(Z_stat)
}

# t1=Sys.time()
# sim_Z_beta()
# t2=Sys.time()
# t2-t1

#################################################################### 


N=1000
ZR<-vector(length=N)


#t1<-Sys.time()

for(t in 1:N)
{
  
  ZR[t]=sim_Z_ER(n = 100, p=0.5, m = 32, r = 4) 
  print(t)
}

#t2<-Sys.time()




tikz("L4RandomGraph.tex", standAlone=TRUE,width=10, height=7)
hist(ZR,freq=F,breaks=20, ylim=c(0, 0.4), cex.axis=1.5, cex.lab=1.5, main="Histogram of $L_4$ Statistic in the ER Model", xlab='$\\hat{Z}_{m, n, 2}$')
lines(seq(-4,4,by=.01),dnorm(seq(-4,4,by=.01),0,1),col="blue", lwd=2)

#qqnorm(Z_vec_pb,cex=0.3,pch=20,col="blue",main=NULL,xlim=c(-4,4),ylim=c(-4,4))
#qqline(Z_vec_pb,col="red")

dev.off()

tools::texi2dvi("L4RandomGraph.tex",pdf=T)


##################################################################### 


#set.seed(2022)
N=1000
ZSBM<-vector(length=N)
  

#t1<-Sys.time()

for(t in 1:N)
{
  
    ZSBM[t]=sim_Z_PB(n = 100, a=8, b=2, m = 32, r = 4) 
    print(t)
}

#t2<-Sys.time()




tikz("L4SBMSparse.tex", standAlone=TRUE,width=10, height=7)
hist(ZSBM,freq=F,breaks=20, ylim=c(0, 0.4), cex.axis=1.5, cex.lab=1.5, main="Histogram of $L_4$ Statistic in the 2-Block SBM", xlab='$\\hat{Z}_{m, n, 4}$')
lines(seq(-4,4,by=.01),dnorm(seq(-4,4,by=.01),0,1),col="blue", lwd=2)

#qqnorm(Z_vec_pb,cex=0.3,pch=20,col="blue",main=NULL,xlim=c(-4,4),ylim=c(-4,4))
#qqline(Z_vec_pb,col="red")

dev.off()

tools::texi2dvi("L4SBMSparse.tex",pdf=T)


##################################################################### 



#set.seed(2022)
N=1000
Zbeta=matrix(0, N, 2)

#t1<-Sys.time()

for(t in 1:N)
{ 
    Zbeta[t,1]=sim_Z_beta(n = 100, rad = 1, m = 32, r = 2)
    Zbeta[t,2]=sim_Z_beta(n = 100, rad = 1, m = 32, r = 4) 
    print(t)
}






tikz("L2Beta.tex", standAlone=TRUE,width=10, height=7)
hist(Zbeta[,1],freq=F,breaks=20, ylim=c(0, 0.4), cex.axis=1.5, cex.lab=1.5, main="Histogram of $L_2$ Statistic in the $\\beta$-Model", xlab='$\\hat{Z}_{m, n, 2}$')
lines(seq(-4,4,by=.01),dnorm(seq(-4,4,by=.01),0,1),col="blue", lwd=2)

#qqnorm(Z_vec_pb,cex=0.3,pch=20,col="blue",main=NULL,xlim=c(-4,4),ylim=c(-4,4))
#qqline(Z_vec_pb,col="red")

dev.off()

tools::texi2dvi("L2beta.tex",pdf=T)




tikz("L4Beta.tex", standAlone=TRUE,width=10, height=7)
hist(Zbeta[,1],freq=F,breaks=20, ylim=c(0, 0.4), cex.axis=1.5, cex.lab=1.5, main="Histogram of $L_4$ Statistic in the $\\beta$-Model", xlab='$\\hat{Z}_{m, n, 4}$')
lines(seq(-4,4,by=.01),dnorm(seq(-4,4,by=.01),0,1),col="blue", lwd=2)

#qqnorm(Z_vec_pb,cex=0.3,pch=20,col="blue",main=NULL,xlim=c(-4,4),ylim=c(-4,4))
#qqline(Z_vec_pb,col="red")

dev.off()

tools::texi2dvi("L4beta.tex",pdf=T)



