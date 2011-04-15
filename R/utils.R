
topKfeats = function(mgr, K, fn="inds1.ff", batchsize=200,
   feat = c("score", "ind", "geneind"), ffind=1, ginds) {
     intests = mgr@fflist[[ffind]]
     if (feat == "score") op = function(x)sort(x, decreasing=TRUE)[1:K]
     else if (feat == "ind") op = function(x)order(x, decreasing=TRUE,na.last=NA)[1:K]
     else if (feat == "geneind") op = function(x)ginds[
                                        order(x,decreasing=TRUE,na.last=NA)[1:K]]
     else stop("feat not recognized")
     tmp = ffrowapply(
       t(apply(intests[i1:i2,],1,op)),
       X=intests, RETURN=TRUE, RETCOL=K,
       BATCHSIZE=batchsize)
     ff(tmp, filename=fn, dim=c(nrow(intests),K), overwrite=TRUE,
       vmode="short")
       }

dotrans = function(spack, snpchr="chr1", K=20, rhs=~sex,
   geneApply=mclapply, shortfac=10, trprefix="trans", batchsize=200,
   genechr=1:22) {
 cs = getSS(spack, snpchr)
 require(annotation(cs), character.only=TRUE)
 cenv = get(paste(gsub(".db", "", annotation(cs)), "CHR", sep=""))
 chtok = gsub("chr", "", snpchr)
 gspl = AnnotationDbi:::mget(as.character(genechr), revmap(cenv))
 if (length(gspl)<2) stop("need genechr length at least 2")
 allg = featureNames(cs)
 indset = lapply(gspl, function(x) match(x, allg))
 tname = function() gsub(".*/|.*\\\\", "", tempfile())
 curtarg = tname()
 run1 = eqtlTests( cs[probeId(gspl[[1]]),], rhs, targdir=curtarg, runname="run1",
   geneApply=geneApply, shortfac=shortfac )
 baseinds = topKfeats( run1, K, fn=paste(trprefix, "baseinds.ff", sep=""), feat="geneind",
    ginds = indset[[1]] )
 basescores = topKfeats( run1, K, fn=paste(trprefix, "basescores.ff", sep=""), feat="score",
    ginds = indset[[1]] )
 for (i in 2:length(gspl) ) {
     cat(i)
     print(date())
     currun = tname()
     nrun = eqtlTests( cs[probeId(gspl[[i]]),], rhs, targdir=curtarg, runname=currun,
         geneApply=geneApply, shortfac=shortfac )
     incrinds = topKfeats( run1, K, fn=paste(trprefix, "incrinds.ff", sep=""), feat="geneind",
        ginds = indset[[i]] )
     incrscores = topKfeats( run1, K, fn=paste(trprefix, "incrscores.ff", sep=""), feat="score",
        ginds = indset[[i]] )
     updateKfeats( basescores, incrscores, baseinds, incrinds, batchsize=batchsize )
     system("rm -rf incrinds.ff")
     system("rm -rf incrscores.ff")
     system(paste("rm -rf ", filename(nrun@fflist[[1]])))
     }
 list(baseinds=baseinds, basescores=basescores)
}

updateKfeats = function( sco1, sco2, ind1, ind2, batchsize=200 ) {
#
# will overwrite sco1, ind1 with the improvements available in sco2, ind2
   snchunk = chunk(1, nrow(sco1), by=batchsize)
   K=ncol(sco1)
   for (i in 1:length(snchunk)) {
      i1 = snchunk[[i]][1]
      i2 = snchunk[[i]][2]
      scos = cbind(sco1[i1:i2,], sco2[i1:i2,])
      ginds = cbind(ind1[i1:i2,], ind2[i1:i2,])
      rowwiseExtract = function(x,y) t(sapply(1:nrow(x), function(row) x[row,][y[row,]]))
      chind = t(apply(scos, 1, function(x)order(x,decreasing=TRUE,na.last=NA)[1:K]))
      sco1[i1:i2,] = rowwiseExtract( scos, chind )
      ind1[i1:i2,] = rowwiseExtract( ginds, chind )
      }
   invisible(NULL)
}

