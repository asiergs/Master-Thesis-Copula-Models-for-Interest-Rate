library(copula)

set.seed(1000)

# normal copula example ----

cop <- normalCopula(param = c(0.8,-0.7,0.6), dim = 3, dispstr = "un")
u_cop <- rCopula(1219,copula = cop)
varnames = expression('U'[1],'U'[2],'U'[3])
colnames(u_cop) <- varnames

svg(file = 'plots/example_norm_cop_U.svg', width = 5, height = 5)
pairs(u_cop, labels = varnames, cex = 0.5, oma = c(2,2,2,2))
dev.off()

# t-Student copula example ----

cop <- tCopula(param = c(0.8,-0.7,0.6),df = 4, dim = 3, dispstr = "un")
u_cop <- rCopula(1219,copula = cop)

svg(file = 'plots/example_tStud_cop_U.svg', width = 5, height = 5)
pairs(u_cop, labels = varnames, cex = 0.5, oma = c(2,2,2,2))
dev.off()

# AMH copula example ----

cop1 <- amhCopula(param = 0.6)
u_cop1 <- rCopula(1219,copula = cop1)
cop2 <- tCopula(param = 0.75)
u_cop2 <- rCopula(1219,copula = cop2)

svg(file = 'plots/example_AMH1_cop_U.svg', width = 5, height = 5)
pairs(u_cop1, labels = varnames[1:2], cex = 0.5, oma = c(2,2,2,2))
dev.off()
svg(file = 'plots/example_AMH2_cop_U.svg', width = 5, height = 5)
pairs(u_cop2, labels = varnames[1:2], cex = 0.5, oma = c(2,2,2,2))
dev.off()

# Clayton copula example ----

cop1 <- claytonCopula(param = 1)
u_cop1 <- rCopula(1219,copula = cop1)
cop2 <- claytonCopula(param = 5)
u_cop2 <- rCopula(1219,copula = cop2)

svg(file = 'plots/example_Clayton1_cop_U.svg', width = 5, height = 5)
pairs(u_cop1, labels = varnames[1:2], cex = 0.5, oma = c(2,2,2,2))
dev.off()
svg(file = 'plots/example_Clayton2_cop_U.svg', width = 5, height = 5)
pairs(u_cop2, labels = varnames[1:2], cex = 0.5, oma = c(2,2,2,2))
dev.off()

# Frank copula example ----

cop1 <- frankCopula(param = -3)
u_cop1 <- rCopula(1219,copula = cop1)
cop2 <- frankCopula(param = 5)
u_cop2 <- rCopula(1219,copula = cop2)

svg(file = 'plots/example_Frank1_cop_U.svg', width = 5, height = 5)
pairs(u_cop1, labels = varnames[1:2], cex = 0.5, oma = c(2,2,2,2))
dev.off()
svg(file = 'plots/example_Frank2_cop_U.svg', width = 5, height = 5)
pairs(u_cop2, labels = varnames[1:2], cex = 0.5, oma = c(2,2,2,2))
dev.off()

# Gumbel-Hougaard copula example ----

cop1 <- gumbelCopula(param = 3)
u_cop1 <- rCopula(1219,copula = cop1)
cop2 <- gumbelCopula(param = 5)
u_cop2 <- rCopula(1219,copula = cop2)

svg(file = 'plots/example_Gumbel1_cop_U.svg', width = 5, height = 5)
pairs(u_cop1, labels = varnames[1:2], cex = 0.5, oma = c(2,2,2,2))
dev.off()
svg(file = 'plots/example_Gumbel2_cop_U.svg', width = 5, height = 5)
pairs(u_cop2, labels = varnames[1:2], cex = 0.5, oma = c(2,2,2,2))
dev.off()

# Joe copula example ----

cop1 <- joeCopula(param = 2)
u_cop1 <- rCopula(1219,copula = cop1)
cop2 <- joeCopula(param = 5)
u_cop2 <- rCopula(1219,copula = cop2)

svg(file = 'plots/example_Joe1_cop_U.svg', width = 5, height = 5)
pairs(u_cop1, labels = varnames[1:2], cex = 0.5, oma = c(2,2,2,2))
dev.off()
svg(file = 'plots/example_Joe2_cop_U.svg', width = 5, height = 5)
pairs(u_cop2, labels = varnames[1:2], cex = 0.5, oma = c(2,2,2,2))
dev.off()

# Nested copula example ----

theta <- c(4,7)
nacList <- list(theta[1],1,list(list(theta[2],c(2,3))))
nestc <- onacopulaL("Frank",nacList)
u_cop <- rCopula(n = 1250,copula = nestc)

svg(file = 'plots/example_nested_cop_U.svg', width = 5, height = 5)

col1 = "black"
col2 = "blue"
par(bg = col1, col = col2, col.axis = col2, col.lab = col2,
    col.main = col2, col.sub = col2, fg = col2)
pairs(u_cop, labels = varnames, cex = 0.5, oma = c(2,2,2,2))
dev.off()



cop <- normalCopula(param = c(0.8,0.7,0.6), dim = 3, dispstr = "un")
u_cop <- rCopula(1219*30,copula = cop)
varnames = expression('U'[1],'U'[2],'U'[3])
colnames(u_cop) <- varnames

dir.create("gifs")

n <- 10
for (i in 1:1000){
  jpeg(file = paste0('gifs/norm_cop',i,'.jpeg'), width = 500, height = 500)
  pairs(u_cop[(i*n+1):(i*n+1219),], labels = varnames, cex = 1, oma = c(2,2,2,2))
  dev.off()
}
