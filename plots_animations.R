# PC adjust in time ----

line2user <- function(line, side) {
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(0:1, 'inches', 'user'))
  y_off <- diff(grconvertY(0:1, 'inches', 'user'))
  switch(side,
         `1` = par('usr')[3] - line * y_off * lh,
         `2` = par('usr')[1] - line * x_off * lh,
         `3` = par('usr')[4] + line * y_off * lh,
         `4` = par('usr')[2] + line * x_off * lh,
         stop("side must be 1, 2, 3, or 4", call.=FALSE))
}

n <- length(swaps_df4$date)

for (i in 1:n){
  date <- as.character(swaps_df4$date[i])
  png(file = paste0("gifs/PCA2_",i,".png"), width = 1200, height = 1200, units = "px", res = 150)
  par(mfrow=c(2,2), oma = c(2,2,2,2), mar = c(1.5,1.5,4,1.5))
  for (i in 1:4) {
    ap <- recordPlot()
    plot(times,swaps_matrix[date,1:dim], type="p", lty = 1, pch = 19,
         lwd = 2, ylim = c(-1,6.5), ylab = "swap", main = paste0(i," PC"), xlab = "years")
      lines(times,swaps_matrix[date,1:dim], type="l", lty = 1,
          lwd = 2)
    lines(times,swaps_red_list[[i]][date,1:dim], type = "l", lty = 2, lwd = 2, col = "red")
    legend("bottomright",legend = c("Real", "PC approx"), lty = 1:2, lwd = c(2,2), ncol=2,
           col= c("black","red"))
    if (i == 2) text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
                     line2user(line=2, side=3), date, xpd=NA, cex=2, font=2)
  }
  dev.off()
}

# 

par(old.par)
svg(file = "plots/ts_PC.svg", width = 10, height = 3)
par(mfrow=c(1,1), oma = c(1,1,1,1), mar = c(4,1,1,1))
plot(dates, tsPC1, type = "l", col = "red", xlab = "time", ylab = "")
lines(dates, tsPC2, type = "l",
      col = "blue")
lines(dates, tsPC3, type = "l",
      col = "green")
legend("bottomright",legend = paste("PC",1:3), col = c("red", "blue", "green"),
       lwd = rep(1,3),ncol=3)
dev.off()

par(old.par)
svg(file = "plots/ts_DJSTOXX.svg", width = 10, height = 3)
par(mfrow=c(1,1), oma = c(1,1,1,1), mar = c(4,1,1,1))
plot(dates, DJSTOXX, type = "l", col = "blue", xlab = "time", ylab = "")
legend("bottomright",legend = paste("EUROSTOXX 600"), col = c("blue"),
       lwd = rep(1),ncol=1)
dev.off()

par(old.par)
svg(file = "plots/ts_USDEUR.svg", width = 10, height = 3)
par(mfrow=c(1,1), oma = c(1,1,1,1), mar = c(4,1,1,1))
plot(dates, USDEUR, type = "l", col = "blue", xlab = "time", ylab = "")
legend("bottomright",legend = paste("USDEUR"), col = c("blue"),
       lwd = rep(1),ncol=1)
dev.off()

