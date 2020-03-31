source("fastqc.R")

qc.files <- list.files("data", pattern="zip", full.names=TRUE)

qc.data <- lapply(qc.files, read.fastqc )

par(cex.axis=1.2)
par(cex.lab=1.2)
cols <- rep(1, sum(b))
cols[ grepl("185", qc.files) ] <- 2
cols[ grepl("190AD", qc.files) ] <- 3
cols[ grepl("190BD", qc.files) ] <- 4
lty <- ifelse(grepl('1P', qc.files), 1, 3)

layout(matrix(c(1,1, 2,3, 4,4), nrow=3, byrow=TRUE))

plot.pos.quals(qc.data, lwd=2, cols=cols, lty=lty)
legend('bottomleft', legend=c('185BD', '190AD', '190BD'), lty=1, col=c(2,3,4), lwd=2, bty='n')
legend(x=1, y=31, legend=c('read 1', 'read 2'), lty=c(1,3), lwd=2, bty='n')

plot.seq.qual(qc.data, lwd=2, cols=cols, lty=lty)
legend('topleft', legend=c('185BD', '190AD', '190BD'), lty=1, col=c(2,3,4), lwd=2, bty='n')
legend(x=par("usr")[1], y=1.2e6, legend=c('read 1', 'read 2'), lty=c(1,3), lwd=2, bty='n')

plot.length.dist(qc.data, lwd=2, cols=cols, lty=lty)
legend('topleft', legend=c('185BD', '190AD', '190BD'), lty=1, col=c(2,3,4), lwd=2, bty='n')
legend(x=par("usr")[1], y=6e5, legend=c('read 1', 'read 2'), lty=c(1,3), lwd=2, bty='n')

ad.names <- plot.adapter.content(qc.data, lwd=2, lty=lty, cols=cols, cex=0.8)
legend('topleft', legend=ad.names, pch=1:length(ad.names), cex=0.8)
legend(x=par("usr")[1], y=1.3, legend=c('185BD', '190AD', '190BD'), lty=1, col=c(2,3,4), lwd=2, bty='n')
legend(x=par("usr")[1], y=1.0, legend=c('read 1', 'read 2'), lty=c(1,3), lwd=2, bty='n')

