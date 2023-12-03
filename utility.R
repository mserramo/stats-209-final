## puts the results of a pair match in a nice form
## Usage: summarize.match(dataset,pairmatch_output)
summarize.match <- function(dat, ms, ps.name="prop", keep.mset=FALSE) {
  adat <- dat
  adat$mset <- ms
  adat <- adat[!is.na(adat$mset),]
  adat.treat <- adat[adat$z==1, ]
  adat.ctrl <- adat[adat$z==0, ]
  
  adat.m <- merge(adat.treat, adat.ctrl, by="mset", suffixes=c(".1", ".0"))
  
  if(!keep.mset) {
    adat.m <- adat.m[, -which(names(adat.m) %in% c("z.1", "z.0", "mset"))]
  } else {
    adat.m <- adat.m[, -which(names(adat.m) %in% c("z.1", "z.0"))]        
  }
  adat.m <- adat.m[, sort(names(adat.m), index.return=TRUE)$ix]
  
  p0.name <- paste0(ps.name,".", 0)
  p1.name <- paste0(ps.name,".",1)
  
  adat.m.tmp.1 <- adat.m[, -which(names(adat.m) %in% c(p0.name, p1.name))]
  adat.m.tmp.2 <- adat.m[, c(p0.name, p1.name)]
  
  adat.m <- cbind(adat.m.tmp.1, adat.m.tmp.2)
  
  return(adat.m)
}