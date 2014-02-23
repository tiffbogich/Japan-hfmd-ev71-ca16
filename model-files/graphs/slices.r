require(plom)

s <- get.settings()

pdf('slice.pdf')
layout(matrix(1:20, 4,5))
par(mar=c(4,2,1,1))
for(pg in c('r0_1:all','r0_2:all','v:all','q:all','e:all','d:all','z:all','sigma:all','iota_1:all','iota_2:all','SS:all','IS:all','SI:all','SR:all','RS:all','rep:EV71_JAP__IASR__inc','rep:CA16_JAP__IASR__inc')){
  sp <- strsplit(pg, ':')[[1]]
  plot.slice(s, res='slice', par=sp[1], group=sp[2])
}
dev.off()
