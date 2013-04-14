require(grid)
vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange_ggplot2 <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
    dots <- list(...)
    n <- length(dots)
    if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
    if(is.null(nrow)) { nrow = ceiling(n/ncol)}
    if(is.null(ncol)) { ncol = ceiling(n/nrow)}
        ## NOTE see n2mfrow in grDevices for possible alternative
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
    ii.p <- 1
    for(ii.row in seq(1, nrow)){
    ii.table.row <- ii.row
    if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
        for(ii.col in seq(1, ncol)){
            ii.table <- ii.p
            if(ii.p > n) break
            print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
            ii.p <- ii.p + 1
        }
    }
}

library(ggplot2)
library(gridExtra)
atable <- read.delim('alphas.tsv')
atable$a <- with(atable, 1-pN*dS/(pS*dN))
atable <- atable[is.finite(atable$a),]
c0 <- atable[atable$clone_id==0,]
c0dN <- data.frame(gen=c0$generation, count=c0$dN, count_type="dN")
c0pN <- data.frame(gen=c0$generation, count=c0$pN, count_type="pN")
c0dS <- data.frame(gen=c0$generation, count=c0$dS, count_type="dS")
c0pS <- data.frame(gen=c0$generation, count=c0$pS, count_type="pS")
c0counts <- data.frame(rbind(c0dN,c0dS,c0pN,c0pS))
# plt <- ggplot(atable,aes(x=generation,y=a))
# for (clone in unique(atable$clone_id)) {
#     clone_as <- subset(atable,atable$clone_id==clone)
#     plt <- plt + geom_smooth(data=clone_as,color=clone)
# }
aplt <- qplot(data=atable,generation,a,
              colour=as.factor(clone_id),
              geom="smooth",
              main="Alpha vs Generation")
c0plt <- qplot(data=c0counts, gen, count,
              colour=as.factor(count_type),
              geom="smooth",
              main="Site Count vs Generation in Clone 0")
plt <- arrange_ggplot2(aplt,c0plt,ncol=1)
print(plt)
