source("sim_gene/ref/simu.R")

simu_df <- data.frame(ages=rep(NA, 2000),
                      optimise=c(rep("TRUE", 1000), rep("FALSE", 1000)))
for (i in 1:1000) {
  r = simu(n=100, rho = 5, delta = 1, blocks = c(100), optimise = T)
  simu_df$ages[i] <- tail(r$ages, n=1)
  r = simu(n=100, rho = 5, delta = 1, blocks = c(100), optimise = F)
  simu_df$ages[i+1000] <- tail(r$ages, n=1)
}
boxplot(ages ~ optimise, data = simu_df)
mean(simu_df$ages[1:1000])
mean(simu_df$ages[1001:2000])
