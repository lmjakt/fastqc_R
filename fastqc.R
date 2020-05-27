## Functions for getting data out of zipped fastqc output
## and plotting them.

## get data from a fastqc analysis in order to automate analysis of large
## numbers of data points.
## WARNING; gets warning about closing unused connections. I am not
## sure how to close these manually as this
read.fastqc <- function(fname){
    if.name <- sub(".+?/([^/]+$)", "\\1", fname)
    if.name <- paste( sub("\\.zip", "", if.name), "/fastqc_data.txt", sep="")
    con <- unz(fname, if.name)
    open(con)
    data.lines <- readLines( con )
    start.i <- grep("^>>", data.lines)
    end.i <- start.i[ grep("^>>END_MODULE", data.lines[start.i]) ]
    start.i <- setdiff(start.i, end.i)
    sec.names <- sub("^>>", "", data.lines[start.i])
    sec.grades <- sub(".+?\t(.+)$", "\\1", sec.names)
    sec.names <- sub("(.+?)\t.+$", "\\1", sec.names)
    start.i <- start.i + 2
    while( sum( grepl( "^#", data.lines[start.i])) > 0){
        b <- grepl( "^#", data.lines[start.i] )
            start.i[ b ] <- start.i[b] + 1
    }
    if(length(start.i) != length(end.i))
        stop("Unmatched section delimiters")
    ## if( any(start.i >= end.i))
    ##     stop("Unreasonable section positions")
##
    tables <- vector(mode='list', length=length(start.i))
    names(tables) <- sec.names
    close(con)
    for(i in 1:length(start.i)){
        con <- unz(fname, if.name)
        open(con)
        if(end.i[i] > start.i[i]){
            tables[[i]] <- read.table( con, skip=start.i[i] - 1,
                                      nrows=(end.i[i] - start.i[i]),
                                      sep="\t", stringsAsFactors=FALSE )
            header <- sub("^#", "", data.lines[start.i[i]-1])
            colnames(tables[[i]]) <- strsplit( header, split="\t" )[[1]]
        }
        if(isOpen(con))
            close(con)
    }
    c(tables, list('grade'=sec.grades))
}

def.colors <- function(l, m=0.75, s=0.8, v=0.8){ hsv( m * 1:l/l, s, v )}

extract.ranges <- function(labs){
    beg <- as.numeric( sub("^([0-9]+)-[0-9]+", "\\1", labs) )
    end <- as.numeric( sub("^[0-9]+-([0-9]+)", "\\1", labs) )
    mid <- (beg + end) / 2
    i <- 1:length(labs)
    list(lab=labs, beg=beg, end=end, mid=mid, i=i)
}


## fqc is a list of lists.. of tables
plot.pos.quals <- function(fqc,  ppar='Mean',
                           mar=c(7.1, 4.1, 4.1, 2.1),
                           cols=def.colors(length(fqc)),
                           lty=rep(1, length(fqc)), ...){
    tname <- "Per base sequence quality"
    ## note that the actual names will also have ": pass", ": fail"
    ## after this..
    x <- fqc[[1]][[tname]][,1]
    if(!all(sapply(fqc, function(y){
        all(x == y[[tname]][,1]) })))
        stop("Different length ranges")
    y <- sapply(fqc, function(m){ m[[tname]][,ppar] })
    par(mar=mar)
    plot(1:length(x), y[,1], xlab='', ylab='mean quality', ylim=range(y), type='n', xaxt='n', xaxs='i',
         main=tname)
    invisible( sapply(1:ncol(y), function(i){ lines(1:length(x), y[,i], col=cols[i], lty=lty[i], ... ) }))
    axis(1, at=1:length(x), labels=x, las=2)
    invisible(list(length=matrix(x, nrow=nrow(y), ncol=ncol(y)), y=y,
                   x=matrix(1:nrow(y), nrow=nrow(y), ncol=ncol(y)),
                   i=matrix(1:ncol(y), nrow=nrow(y), ncol=ncol(y), byrow=TRUE)))
}

plot.seq.qual <- function(fqc, mar=c(5.1, 4.1, 4.1, 2.1),
                          cols=def.colors(length(fqc)),
                          lty=rep(1, length(fqc)), ...){
    tname <- "Per sequence quality scores"
    x <- lapply(fqc, function(m){ m[[tname]][,'Quality'] })
    y <- lapply(fqc, function(m){ m[[tname]][,'Count'] })
    plot(1, 1, type='n', xlab='Mean base quality', ylab='Count',
         main=tname, xlim=range(unlist(x)), ylim=range(unlist(y)) )
    invisible(sapply(1:length(x), function(i){
        lines( x[[i]], y[[i]], col=cols[i], lty=lty[i], ... )}))
}

plot.length.dist <- function(fqc,
                             mar=c(7.1, 4.1, 4.1, 2.1),
                             cols=def.colors( length(fqc) ),
                             lty=rep(1, length(fqc)), ...){
    tname <- "Sequence Length Distribution"
    ## note that the actual names will also have ": pass", ": fail"
    ## after this..
    x <- lapply( fqc, function(m){
        labs <- m[[tname]][,'Length']
        beg <- as.numeric( sub("^([0-9]+)-[0-9]+", "\\1", labs) )
        end <- as.numeric( sub("^[0-9]+-([0-9]+)", "\\1", labs) )
        mids <- (beg + end) / 2
        i <- 1:length(labs)
        list(lab=labs, beg=beg, end=end, mid=mids, i=i)
    })
    y <- lapply(fqc, function(m){ m[[tname]][,'Count'] })    
    par(mar=mar)
    plot(unlist(sapply(x, function(z){ z$mid })),
         unlist(y), xlab='', ylab='Count', type='n', xaxt='n', xaxs='i',
         main=tname)
    invisible( lapply(1:length(fqc),
                      function(i){ lines(x[[i]]$mid,
                                         y[[i]],
                                         col=cols[i], lty=lty[i], ... ) }))
    x.x <- unlist( lapply( x, function(z){ z$mid }))
    x.lab <- unlist( lapply( x, function(z){ z$lab }))
    x.lab <- x.lab[ !duplicated(x.x) ]
    x.x <- x.x[ !duplicated(x.x) ]
    axis(1, at=x.x, labels=x.lab, las=2)
}

plot.adapter.content <- function(fqc,
                             mar=c(7.1, 4.1, 4.1, 2.1),
                             cols=def.colors(length(fqc)),
                             lty=rep(1, length(fqc)), ylim=NULL,...){
    tname <- "Adapter Content"
    ## note that the actual names will also have ": pass", ": fail"
    ## after this..
    x <- fqc[[1]][[tname]][,1]
    if(!all(sapply(fqc, function(y){
        all(x == y[[tname]][,1]) })))
        stop("Different length ranges")
    y <- lapply(colnames(fqc[[1]][[tname]])[-1], function(column){
        sapply(fqc, function(m){ m[[tname]][,column] }) })
    par(mar=mar)
    if(is.null(ylim))
        ylim <- range( unlist(y) )
    plot(1:length(x), y[[1]][,1], xlab='', ylab='percent', ylim=ylim, type='n', xaxt='n', xaxs='i',
         main=tname)
    invisible( lapply(1:length(y), function(i){
        yy <- y[[i]]
        sapply(1:ncol(yy), function(j){ lines(1:length(x), yy[,j], col=cols[j], lty=lty[j], type='l', pch=i, ... ) }) }))
    invisible( lapply(1:length(y), function(i){
        yy <- y[[i]]
        sapply(1:ncol(yy), function(j){ lines(1:length(x), yy[,j], col=cols[j], lty=lty[j], type='p', pch=i, ... ) }) }))
    axis(1, at=1:length(x), labels=x, las=2)
    invisible(colnames(fqc[[1]][[tname]])[-1])
}

plot.gc <- function(fqc,
##                    mar=c(7.1, 4.1, 4.1, 2.1),
                    cols=def.colors( length(fqc) ),
                    lty=rep(1, length(fqc)), ...){
    gc.pcnt <- lapply(fqc, function(x){ x[['Per sequence GC content']][,1] })
    gc.count <- lapply(fqc, function(x){ x[['Per sequence GC content']][,2] })
    plot(1,1, type='n', xlab='GC%', ylab='count',
         xlim=c(0,100), ylim=range( unlist(gc.count)), ... )
    for(i in 1:length(fqc)){
        lines( gc.pcnt[[i]], gc.count[[i]], col=cols[i], lty=lty[i], ...)
    }
    invisible(list(pc=gc.pcnt, count=gc.count))
}

plot.bases <- function(fqc,
                       cols=def.colors( length(fqc) ),
                       lty=rep(1, length(fqc)), ...){
    tbl <- lapply( fqc, function(x){ x[['Per base sequence content']] })
    x <- extract.ranges( unlist( lapply( tbl, function(y){ y[,'Base'] }) ))
    pcnt <- unlist( lapply( tbl, function(y){ y[,-1] }) )
    for(i in 2:5){
        plot(1,1, type='n', xlab='position', ylab='base %',
             xlim=range(x$mid), ylim=c(0, max(pcnt)) )
        for(j in 1:length(tbl)){
            xx <- extract.ranges( tbl[[j]][,1] )
            lines( xx$mid, tbl[[j]][,i], col=cols[j], lty=lty[j] )
        }
    }
}

plot.dup <- function(fqc,
                     cols=def.colors( length(fqc) ),
                     lty=rep(1, length(fqc)), ...){
    dup <- lapply(fqc, function(x){ x[['Sequence Duplication Levels']] })
    dup.l <- unlist( lapply(dup, function(x){ x[,1] }) )
    pcn.1 <- unlist( lapply(dup, function(x){ x[,2] }) )
    pcn.2 <- unlist( lapply(dup, function(x){ x[,3] }) )
    plot(1,1, type='n', xlab='Duplication level', ylab='% of deduplicated',
         xlim=c(1,nrow(dup[[1]])), ylim=c(0,max(pcn.1)), axes=FALSE )
    for(i in 1:length(dup))
        lines( 1:nrow(dup[[i]]), dup[[i]][,2], col=cols[i], lty=lty[i], ... )
    axis(1, at=1:nrow(dup[[1]]), labels=dup[[1]][,1])
    axis(2)
    
    plot(1,1, type='n', xlab='Duplication level', ylab='% of total',
         xlim=c(1,nrow(dup[[1]])), ylim=c(0,max(pcn.2)), axes=FALSE )
    axis(1, at=1:nrow(dup[[1]]), labels=dup[[1]][,1])
    axis(2)
    for(i in 1:length(dup))
        lines( 1:nrow(dup[[1]]), dup[[i]][,3], col=cols[i], lty=lty[i], ... )
}

plot.grades <- function(fqc){
    grades <- sapply( fqc, function(x){ x[['grade']] })
    y <- matrix(1:nrow(grades), nrow=nrow(grades), ncol=ncol(grades))
    x <- matrix(1:ncol(grades), nrow=nrow(grades), ncol=ncol(grades), byrow=T)
    cols <- c( rgb(0, 0.8, 0), rgb(0.5, 0.5, 0), rgb(0.8, 0, 0) )
    names(cols) <- c('pass', 'warn', 'fail')
    plot.new()
    plot.window(xlim=c(1, max(x)+1), ylim=c(1, max(y)+1))
    rect( x, y, x+1, y+1, col=cols[grades] )
}
       
plot.read.n <- function(fqc, sm.lab=NULL, ...){
    counts <- sapply(fqc, function(x){
        m <- x[['Basic Statistics']]
        i <- m[,1] == 'Total Sequences'
        as.numeric(m[i, 'Value'])})
    barplot( counts, names.arg=sm.lab, ... )
    invisible( counts )
}





