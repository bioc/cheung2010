library(multicore)
library(cheung2010)

doalltrans = function(spack, snpchr="chr20", K=25, rhs=~sex,
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
 outmgr = list()
 indslist = list()
 scoreslist = list()
 for (i in 1:length(gspl)) {
   cat(i)
 
   outmgr[[i]] = eqtlTests( cs[probeId(gspl[[i]]),], rhs, targdir=curtarg, runname=paste("gc", i, sep=""), geneApply=geneApply, shortfac=shortfac )
   indslist[[i]] = topKfeats( outmgr[[i]], K=K, fn=paste("inds",trinfix, i, ".ff", sep=""), feat="geneind",
     ginds = indset[[i]] )
   scoreslist[[i]] = topKfeats( outmgr[[i]], K=K, fn=paste("scores", trinfix, i, ".ff", sep=""), feat="score",
     ginds = indset[[i]] )
 }
 save(outmgr, file="outmgr.rda")
 list(outmgr=outmgr, indslist=indslist, scoreslist=scoreslist, geneuniverse=allg, K=K, trinfix=trinfix)
}

sss = doalltrans("cheung2010", genechr=1:22, trinfix="D")
save(sss, file="sss.rda")

for (j in 2:length(sss$scoreslist) )
 updateKfeats( sss$scoreslist[[1]], sss$indslist[[1]], sss$scoreslist[[j]], sss$indslist[[j]] )
