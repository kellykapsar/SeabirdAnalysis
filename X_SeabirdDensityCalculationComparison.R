

# See notes in working methods onenote for 13 July 2023 for explanation 
nsamples <- 10000

den1list <- c()
den2list <- c()

for(i in 1:nsamples){
  sobs <- data.frame(sample_area = runif(1000, 0.7, 0.9),
                     nbirds = rnbinom(nsamples, 10, .5))
  
  den1 <- sum(sobs$nbirds)/sum(sobs$sample_area) # total estimates 
  den2 <- mean(sobs$nbirds/sobs$sample_area) # Transect-based estimates 
  den1list <- c(den1list, den1)
  den2list <- c(den2list, den2)
}

summary(den1list)
summary(den2list)

df <- data.frame(method=as.factor(c(rep("tot_avg", nsamples), rep("tran_avg", nsamples))), estimate=c(den1list, den2list))

boxplot(df$estimate ~ df$method,
        col='steelblue',
        main='Density estimate comparison',
        xlab='Method',
        ylab='Density estimate') 
                 