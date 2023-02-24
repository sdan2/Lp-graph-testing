

powerrandomgraph<-read.table(file="RandomGraphSeparationRevised.txt")


pdf(file="RandomGraphSeparationRevised.pdf")

plot(powerrandomgraph[,1], powerrandomgraph[,3], type='b', col=1, pch=1, lwd=1.5, ylim=c(0, 1), xlab="Separation", ylab="Power", main="Power in the ER Model for Sparse Alternatives")
points(powerrandomgraph[,1], powerrandomgraph[,4], type='b', col=2, pch=2, lwd=1.5) 
points(powerrandomgraph[,1], powerrandomgraph[,5], type='b', col=3, pch=3, lwd=1.5) 
points(powerrandomgraph[,1], powerrandomgraph[,6], type='b', col=4, pch=4, lwd=1.5) 


legend("bottomright", c("L2 Norm", "L4 Norm", "L6 Norm", "Max Norm"), col=c(1, 2, 3, 4), pch = c(1, 2, 3, 4), bg = 'gray90')

dev.off() 




###############################################


powerrandomgraph<-read.table(file="SBM 10 Separation.txt")


pdf(file="SBM10Separation.pdf")

plot(powerrandomgraph[,1], powerrandomgraph[,3], type='b', col=1, pch=1, lwd=1.5, ylim=c(0, 1), xlab="Separation", ylab="Power", main="Power in the PBM for Sparse Alternatives")
points(powerrandomgraph[,1], powerrandomgraph[,4], type='b', col=2, pch=2, lwd=1.5) 
points(powerrandomgraph[,1], powerrandomgraph[,5], type='b', col=3, pch=3, lwd=1.5) 


legend("bottomright", c("L2 Norm", "L4 Norm", "Max Norm"), col=c(1, 2, 3), pch = c(1, 2, 3), bg = 'gray90')

dev.off() 


###################################################


powerrandomgraph<-read.table(file="SBM 50 Separation.txt")


pdf(file="SBM50Separation.pdf")

plot(powerrandomgraph[,1], powerrandomgraph[,3], type='b', col=1, pch=1, lwd=1.5, ylim=c(0, 1), xlab="Separation", ylab="Power", main="Power in the PBM for Sparse Alternatives")
points(powerrandomgraph[,1], powerrandomgraph[,4], type='b', col=2, pch=2, lwd=1.5) 
points(powerrandomgraph[,1], powerrandomgraph[,5], type='b', col=3, pch=3, lwd=1.5) 


legend("bottomright", c("L2 Norm", "L4 Norm", "Max Norm"), col=c(1, 2, 3), pch = c(1, 2, 3), bg = 'gray90')

dev.off() 



###############################################

