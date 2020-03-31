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
    if( any(start.i >= end.i))
        stop("Unreasonable section positions")
##
    tables <- vector(mode='list', length=length(start.i))
    names(tables) <- sec.names
    close(con)
    for(i in 1:length(start.i)){
        con <- unz(fname, if.name)
        tables[[i]] <- read.table( con, skip=start.i[i] - 1,
                          nrows=(end.i[i] - start.i[i]), sep="\t", stringsAsFactors=FALSE )
        header <- sub("^#", "", data.lines[start.i[i]-1])
        colnames(tables[[i]]) <- strsplit( header, split="\t" )[[1]]
    }
    c(tables, list('grade'=sec.grades))
}

def.colors <- function(l, m=0.75, s=0.8, v=0.8){ hsv( m * 1:l/l, s, v )}

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
    x <- sapply( fqc, function(m){ m[[tname]][,'Length'] })
    y <- sapply(fqc, function(m){ m[[tname]][,'Count'] })    
    par(mar=mar)
    plot(x, y, xlab='', ylab='Count', ylim=range(y), type='n', xaxt='n', xaxs='i',
         main=tname)
    invisible( sapply(1:length(fqc),
                      function(i){ lines(fqc[[i]][[tname]][,'Length'],
                                         fqc[[i]][[tname]][,'Count'],
                                         col=cols[i], lty=lty[i], type='b', ... ) }))
    x <- unique(x)
    axis(1, at=1:length(x), labels=x, las=2)
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






