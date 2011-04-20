
topKfeats = function(mgr, K, fn="inds1.ff", batchsize=200,
   feat = c("score", "ind", "geneind"), ffind=1, ginds) {
     intests = mgr@fflist[[ffind]]
     if (feat == "score") op = function(x)sort(x, decreasing=TRUE)[1:K]
     else if (feat == "ind") op = function(x)order(x, decreasing=TRUE)[1:K]
     else if (feat == "geneind") op = function(x)ginds[
                                        order(x,decreasing=TRUE)[1:K]]
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
     system(paste("rm -rf ", trprefix, "incrinds.ff", sep=""))
     system(paste("rm -rf ", trprefix, "incrscores.ff", sep=""))
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
      chind = t(apply(scos, 1, function(x)order(x,decreasing=TRUE)[1:K]))
      sco1[i1:i2,] = rowwiseExtract( scos, chind )
      ind1[i1:i2,] = rowwiseExtract( ginds, chind )
      }
   invisible(NULL)
}

transScores = function(smpack, snpchr="chr1", rhs, K=20, targdir="tsco", geneApply=mclapply,
   chrnames=as.character(1:22)) {
 if (length(chrnames)<2) stop("must have length(chrnames) >= 2")
 theCall = match.call()
 require(GGtools)
 sms = getSS(smpack, snpchr)
 guniv = featureNames(sms)
 smanno = gsub(".db", "", annotation(sms))
 require(paste(smanno, ".db", sep=""), character.only=TRUE)
 pnameList = mget(chrnames, revmap(get(paste(smanno, "CHR", sep=""))))
 genemap = lapply( pnameList, function(x) match(x, guniv) )
 nchr = length(chrnames)
 inimgr = eqtlTests( sms[ probeId(pnameList[[chrnames[1]]]),], rhs,
      targdir=targdir, runname=paste("tsc_", chrnames[1], sep=""), geneApply=geneApply,
      saveSummaries=FALSE )
 topKinds = topKfeats( inimgr, K=K, fn=paste(targdir, "/tsinds1_1.ff", sep=""), feat="geneind", ginds=genemap[[1]] )
 topKscores = topKfeats( inimgr, K=K, fn=paste(targdir, "/tssco1_1.ff", sep=""), feat="score", ginds=genemap[[1]] )
 unlink(filename(inimgr@fflist[[1]]))
 for (j in 2:nchr) {  # get scores for same set of SNPs against a new set of genes (next chrom of genes)
    cat(j)
    nxtmgr = eqtlTests( sms[ probeId(pnameList[[chrnames[j]]]),], rhs,
         targdir=targdir, runname=paste("tsctmp", j, sep=""), geneApply=geneApply, saveSummaries=FALSE )
    nxtKinds = topKfeats( nxtmgr, K=K, fn=paste(targdir, "indscratch.ff", sep=""), feat="geneind", ginds=genemap[[j]] )
    nxtKscores = topKfeats( nxtmgr, K=K, fn=paste(targdir, "scoscratch.ff", sep=""), feat="score", ginds=genemap[[j]] )
    updateKfeats( topKscores, nxtKscores, topKinds, nxtKinds )
    unlink(filename(nxtmgr@fflist[[1]]))
    unlink(paste(targdir, "indscratch.ff", sep=""))
    unlink(paste(targdir, "scoscratch.ff", sep=""))
    }
 list(scores=topKscores, inds=topKinds, guniv=guniv, snpnames=rownames(inimgr@fflist[[1]]),
   call=theCall)
}

