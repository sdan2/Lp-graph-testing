library(igraph)
library(Matrix)
library(tikzDevice) 



sim_Z_ER <- function(n = 100, p=0.2, m = 32, r = 4, epsilon)
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
      B = as_adj(sample_gnp(n, min(p+epsilon/n^(2/r), 1))) 
      Diff = Diff + as(Matrix(A - B, sparse = F), "sparseMatrix")
      Sum = Sum + Matrix(A , sparse = F)
    }
    if(l%%2 == 0)
    {
      Prod_stat = Prod_stat*(as(Matrix(Temp + Diff, sparse = F), "sparseMatrix"))
    }
    if(l <= r)
    {
      Prod_var = Prod_var*(Sum/c)
    }
    else
    {
      Prod_var = Prod_var*(1 - Sum/c)
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

sim_Z_PB <- function(n = 100, a=8, b=2, m = 32, r = 4, epsilon)
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
      s1=min(max(a/n^(2/r), 0), 1) 
      s2=min(max(b/n^(2/r), 0), 1)
      t1=min(max(a/n^(2/r)+epsilon/n^(2/r), 0), 1) 
      t2=min(max(b/n^(2/r)-epsilon/n^(2/r), 0), 1)
      A = as_adj(sample_sbm(n, pref.matrix = matrix(c(s1, s2, s2, s1), 2, 2), block.sizes = c(n/2, n/2)))
      B = as_adj(sample_sbm(n, pref.matrix = matrix(c(t1, t2, t2, t1), 2, 2), block.sizes = c(n/2, n/2)))
      Diff = Diff + as(Matrix(A - B, sparse = F), "sparseMatrix")
      Sum = Sum + Matrix(A, sparse = F)
    }
    if(l%%2 == 0)
    {
      Prod_stat = Prod_stat*(as(Matrix(Temp + Diff, sparse = F), "sparseMatrix"))
    }
    if(l <= r)
    {
      Prod_var = Prod_var*(Sum/c)
    }
    else
    {
      Prod_var = Prod_var*(1 - Sum/c)
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

sample_beta <- function(n = 100, beta)
{
  
  #beta = rnorm(n)
  #beta = rad*c(1: n) 
  adj = matrix(0, n, n)
  for(.i in 2:n)
  {                         
    for(.j in 1:(.i-1))
    {
      .p = graphon(beta[.i], beta[.j])
      adj[.i, .j] = rbinom(1, 1, .p)
    }
  }
  adjsymm = adj + t(adj) 
  return(adjsymm)
}

sim_Z_beta <- function(n = 100, m = 32, r = 4, beta, epsilon)
{
  c = m/(2*r)
  b = n/2
  Prod_var = matrix(1, n, n) 
  Prod_stat = matrix(1, n, n) 
  Diff =  matrix(0, n, n) 
  for(l in 1:(2*r))
  {
    Temp = Diff
    Diff =  matrix(0, n, n) 
    Sum =   matrix(0, n, n) 
    for(s in 1:c)
    {
      
      A = sample_beta(n, beta) 
      B = sample_beta(n, beta+rep(epsilon/n^(2/r), n)) 
      Diff = Diff + Matrix(A - B, sparse = F) 
      Sum = Sum + Matrix(A, sparse = F)
    }
    if(l%%2 == 0)
    {
      Prod_stat = Prod_stat*(Matrix(Temp + Diff, sparse = F))
    }
    if(l <= r)
    {
      Prod_var = Prod_var*(Sum/c)
    }
    else
    {
      Prod_var = Prod_var*(1 - Sum/c)
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


#set.seed(2022)
epsilon=seq(0, 5, length.out=10)
ZR<-vector(length=10) 
power2<-vector(length=10)
power4<-vector(length=10)
power6<-vector(length=10)


#t1<-Sys.time()

iterations=100

for(t in 1:10)
{
  
  for(j in 1:iterations)
  {
    ZR[j]=sim_Z_ER(n = 50, p=0.2, m = 48, r = 2, epsilon=epsilon[t]) 
  }
  
  power2[t]<-length(which(abs(ZR)>1.96))/iterations
  
  print(power2[t]) 
}

#t2<-Sys.time()


for(t in 1:10)
{
  
  for(j in 1:iterations)
  {
    ZR[j]=sim_Z_ER(n = 50, p=0.2, m = 48, r = 4, epsilon=epsilon[t]) 
  }
  
  power4[t]<-length(which(abs(ZR)>1.96))/iterations
  
  print(power4[t]) 
  
}



for(t in 1:10)
{
  
  for(j in 1:iterations)
  {
    ZR[j]=sim_Z_ER(n = 50, p=0.2, m = 48, r = 6, epsilon=epsilon[t]) 
  }
  
  power6[t]<-length(which(abs(ZR)>1.96))/iterations
  
  print(power6[t]) 
  
}




power=cbind(power2, power4, power6)
write.table(power, file="RandomGraphSampleSize.txt") 

powerrandomgraph<-read.table(file="RandomGraphSampleSize.txt") 
epsilon=seq(0, 5, length.out=10)


pdf(file="RandomGraphSampleSize.pdf")

plot(epsilon, powerrandomgraph[,1], type='b', col=1, pch=1, lwd=1.5, ylim=c(0, 1), xlab="Separation", ylab="Power", main="Power in the ER Model")
points(epsilon, powerrandomgraph[,2], type='b', col=2, pch=2, lwd=1.5) 
points(epsilon, powerrandomgraph[,3], type='b', col=3, pch=3, lwd=1.5) 


legend("bottomright", c("L2 Norm", "L4 Norm", "L6 Norm"), col=c(1, 2, 3), pch = c(1, 2, 3), bg = 'gray90')

dev.off() 

#################################################################### 



#set.seed(2022) 
epsilon=seq(0, 4, length.out=10)
ZSBM<-vector(length=10) 
power2<-vector(length=10)
power4<-vector(length=10)
power6<-vector(length=10) 

#t1<-Sys.time()

iterations=100


for(t in 1:10)
{
  
  for(j in 1:iterations)
  {
    ZSBM[j]=sim_Z_PB(n = 50, a=2, b=1, m = 48, r = 2, epsilon=epsilon[t]) 
  }
  
  power2[t]<-length(which(abs(ZSBM)>1.96))/iterations
  
  print(power2[t]) 
}

#t2<-Sys.time()


iterations=100 


for(t in 1:10)
{
  
  for(j in 1:iterations)
  {
    ZSBM[j]=sim_Z_PB(n = 50, a=2, b=1, m = 48, r = 4, epsilon=epsilon[t]) 
  }
  
  power4[t]<-length(which(abs(ZSBM)>1.96))/iterations
  
  print(power4[t]) 
}

#t2<-Sys.time()


iterations=100


for(t in 1:10)
{
  
  for(j in 1:iterations)
  {
    ZSBM[j]=sim_Z_PB(n = 50, a=2, b=1, m = 48, r = 6, epsilon=epsilon[t]) 
  }
  
  power6[t]<-length(which(abs(ZSBM)>1.96))/iterations
  
  print(power6[t]) 
}



power=cbind(power2, power4, power6)
write.table(power, file="SBMSampleSize.txt") 


powerSBM<-read.table(file="SBMSampleSize.txt") 
epsilon=seq(0, 4, length.out=10)


pdf(file="SBMSampleSize.pdf")

plot(epsilon, powerSBM[,1], type='b', col=1, pch=1, lwd=1.5, ylim=c(0, 1), xlab="Separation", ylab="Power", main="Power in the PBM")
points(epsilon, powerSBM[,2], type='b', col=2, pch=2, lwd=1.5) 
points(epsilon, powerSBM[,3], type='b', col=3, pch=3, lwd=1.5) 


legend("bottomright", c("L2 Norm", "L4 Norm", "L6 Norm"), col=c(1, 2, 3), pch = c(1, 2, 3), bg = 'gray90')

dev.off() 


#################################################################### 



  
  #set.seed(2022) 
  epsilon=seq(0, 10, length.out=10)
  Zbeta<-vector(length=10) 
  power2<-vector(length=10)
  power4<-vector(length=10)
  power6<-vector(length=10) 
  
  #t1<-Sys.time()
  
  iterations=100
  
  beta=rnorm(n)
  
  for(t in 1:10)
  {
    
    for(j in 1:iterations)
    {
      Zbeta[j]=sim_Z_beta(n = 50, m = 48, r = 2, beta, epsilon=epsilon[t]) 
    }
    
    power2[t]<-length(which(abs(Zbeta)>1.96))/iterations
    
    print(power2[t]) 
  }
  
  #t2<-Sys.time()
  
  


iterations=100


for(t in 1:10)
{
  
  for(j in 1:iterations)
  {
    Zbeta[j]=sim_Z_beta(n = 50, m = 48, r = 4, beta, epsilon=epsilon[t]) 
  }
  
  power4[t]<-length(which(abs(Zbeta)>1.96))/iterations
  
  print(power4[t]) 
}

#t2<-Sys.time()




iterations=100


for(t in 1:10)
{
  
  for(j in 1:iterations)
  {
    Zbeta[j]=sim_Z_beta(n = 50, m = 48, r = 6, beta, epsilon=epsilon[t]) 
  }
  
  power6[t]<-length(which(abs(Zbeta)>1.96))/iterations
  
  print(power6[t]) 
}





power=cbind(power2, power4, power6)
write.table(power, file="BetaModelSampleSize.txt") 


powerbetamodel<-read.table(file="BetaModelSampleSize.txt") 
epsilon=seq(0, 10, length.out=10)


pdf(file="BetaModelSampleSize.pdf")

plot(epsilon, powerbetamodel[,1], type='b', col=1, pch=1, lwd=1.5, ylim=c(0, 1), xlab="Separation", ylab="Power", main="Power in the Beta Model")
points(epsilon, powerbetamodel[,2], type='b', col=2, pch=2, lwd=1.5) 
points(epsilon, powerbetamodel[,3], type='b', col=3, pch=3, lwd=1.5) 


legend("bottomright", c("L2 Norm", "L4 Norm", "L6 Norm"), col=c(1, 2, 3), pch = c(1, 2, 3), bg = 'gray90')

dev.off() 




