
topKfeats = function(mgr, K, fn="inds1.ff", batchsize=200,
   feat = c("score", "ind", "geneind"), ffind=1, ginds ) {
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

.transScores = function(smpack, snpchr="chr1", rhs, K=20, targdir="tsco", geneApply=mclapply,
   chrnames=paste("chr", as.character(1:22), sep=""), geneRanges=NULL, snpRanges=NULL, radius=2e6) {
 if (length(chrnames)<2) stop("must have length(chrnames) >= 2")
 theCall = match.call()
 require(GGtools)
 sms = getSS(smpack, snpchr)
 guniv = featureNames(sms)
 smanno = gsub(".db", "", annotation(sms))
 require(paste(smanno, ".db", sep=""), character.only=TRUE)
 clcnames = gsub("chr", "", chrnames)
 pnameList = mget(clcnames, revmap(get(paste(smanno, "CHR", sep=""))))
 names(pnameList) = chrnames
 genemap = lapply( pnameList, function(x) match(x, guniv) )
 nchr = length(chrnames)
 inimgr = eqtlTests( sms[ probeId(pnameList[[chrnames[1]]]),], rhs,
      targdir=targdir, runname=paste("tsc_", chrnames[1], sep=""), geneApply=geneApply,
      saveSummaries=FALSE )
 if (snpchr == chrnames[1]) {   # kill off scores for genes within radius of SNPs
        if (is.null(geneRanges) || is.null(snpRanges)) stop("ranges must be supplied to exclude cis tests")
        dimnt = dimnames(inimgr@fflist[[1]])
        if (!isTRUE(all(dimnt[[2]] %in% names(geneRanges)))) stop("geneRanges must include names/ranges for all colnames of mgr fflist")
        if (!isTRUE(all(dimnt[[1]] %in% names(snpRanges)))) stop("snpRanges must include names/ranges for all rownames of mgr fflist")
        ol = findOverlaps( snpRanges[ dimnt[[1]] ], geneRanges[ dimnt[[2]] ]+radius )
        matm = matchMatrix( ol )
        if (nrow(matm) > 0) {
             inimgr@fflist[[1]][ matm ] = 0
             } 
        }
 topKinds = topKfeats( inimgr, K=K, fn=paste(targdir, "/", snpchr, "_tsinds1_1.ff", sep=""), feat="geneind", ginds=genemap[[1]] )
 topKscores = topKfeats( inimgr, K=K, fn=paste(targdir, "/", snpchr, "_tssco1_1.ff", sep=""), feat="score", ginds=genemap[[1]] )
 unlink(filename(inimgr@fflist[[1]]))
 for (j in 2:nchr) {  # get scores for same set of SNPs against a new set of genes (next chrom of genes)
    cat(j)
    nxtmgr = eqtlTests( sms[ probeId(pnameList[[chrnames[j]]]),], rhs,
         targdir=targdir, runname=paste("tsctmp", j, sep=""), geneApply=geneApply, saveSummaries=FALSE )
    if (snpchr == chrnames[j]) {   # kill off scores for genes within radius of SNPs
        if (is.null(geneRanges) || is.null(snpRanges)) stop("ranges must be supplied to exclude cis tests")
        dimnt = dimnames(nxtmgr@fflist[[1]])
        if (!isTRUE(all(dimnt[[2]] %in% names(geneRanges)))) stop("geneRanges must include names/ranges for all colnames of mgr fflist")
        if (!isTRUE(all(dimnt[[1]] %in% names(snpRanges)))) stop("snpRanges must include names/ranges for all rownames of mgr fflist")
        ol = findOverlaps( snpRanges[ dimnt[[1]] ], geneRanges[ dimnt[[2]] ]+radius )
        matm = matchMatrix( ol )
        if (nrow(matm) > 0) {
             inimgr@fflist[[1]][ matm ] = 0
             } 
        }
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

cisZero = function(mgr, snpRanges, geneRanges, radius) {
        dimnt = dimnames(mgr@fflist[[1]])
        oks = gsub("rs", "", dimnt[[1]])
        okg = dimnt[[2]]
        srids = values(snpRanges)$RefSNP_id
        if (!isTRUE(all(okg %in% names(geneRanges))))  {
            warning("geneRanges does not include names/ranges for all colnames of mgr fflist")
            okg = intersect( okg, names(geneRanges))
            }
        if (!isTRUE(all(oks %in% srids)))  {
            warning("snpRanges does not include names/ranges for all rownames of mgr fflist")
 # match and numeric indexing much more efficient than names for large objects
            oks = match(oks, srids, nomatch=0)
            if (length(oks) == 0) stop("snpRanges does not include any snps in mgr")
            }
 # can't use zero indices in subsetting GRanges...
        ol = findOverlaps(snpRanges[oks[oks>0]], geneRanges[okg] + 
            radius)
        matm = matchMatrix(ol)
        if (nrow(matm) > 0) {
            mgr@fflist[[1]][matm] = 0
        }
    }

transScores = function (smpack, snpchr = "chr1", rhs, K = 20, targdir = "tsco", 
    geneApply = mclapply, chrnames = paste("chr", as.character(1:22), sep=""), 
    geneRanges = NULL, snpRanges = NULL, radius = 2e+06) 
{
#
# revised with name consistency code
#
    if (length(chrnames) < 2) 
        stop("must have length(chrnames) >= 2")
    theCall = match.call()
    require(GGtools)
    sms = getSS(smpack, snpchr)
    guniv = featureNames(sms)
    smanno = gsub(".db", "", annotation(sms))
    require(paste(smanno, ".db", sep = ""), character.only = TRUE)
    clcnames = gsub("chr", "", chrnames)
    pnameList = mget(clcnames, revmap(get(paste(smanno, "CHR", 
        sep = ""))))
    names(pnameList) = chrnames
    genemap = lapply(pnameList, function(x) match(x, guniv))
    nchr = length(chrnames)
    inimgr = eqtlTests(sms[probeId(pnameList[[chrnames[1]]]), 
        ], rhs, targdir = targdir, runname = paste("tsc_", chrnames[1], 
        sep = ""), geneApply = geneApply, saveSummaries = FALSE)
    if (snpchr == chrnames[1]) {
        if (is.null(geneRanges) || is.null(snpRanges)) 
            stop("ranges must be supplied to exclude cis tests")
        cisZero(inimgr, snpRanges, geneRanges, radius)
    }
    topKinds = topKfeats(inimgr, K = K, fn = paste(targdir, "/", 
        snpchr, "_tsinds1_1.ff", sep = ""), feat = "geneind", 
        ginds = genemap[[1]])
    topKscores = topKfeats(inimgr, K = K, fn = paste(targdir, 
        "/", snpchr, "_tssco1_1.ff", sep = ""), feat = "score", 
        ginds = genemap[[1]])
    unlink(filename(inimgr@fflist[[1]]))
    for (j in 2:nchr) {
        cat(j)
        nxtmgr = eqtlTests(sms[probeId(pnameList[[chrnames[j]]]), 
            ], rhs, targdir = targdir, runname = paste("tsctmp", 
            j, sep = ""), geneApply = geneApply, saveSummaries = FALSE)
        if (snpchr == chrnames[j]) {
            if (is.null(geneRanges) || is.null(snpRanges)) 
                stop("ranges must be supplied to exclude cis tests")
            cisZero(nxtmgr, snpRanges, geneRanges, radius)
            }
        nxtKinds = topKfeats(nxtmgr, K = K, fn = paste(targdir, 
            "indscratch.ff", sep = ""), feat = "geneind", ginds = genemap[[j]])
        nxtKscores = topKfeats(nxtmgr, K = K, fn = paste(targdir, 
            "scoscratch.ff", sep = ""), feat = "score", ginds = genemap[[j]])
        updateKfeats(topKscores, nxtKscores, topKinds, nxtKinds)
        unlink(filename(nxtmgr@fflist[[1]]))
        unlink(paste(targdir, "indscratch.ff", sep = ""))
        unlink(paste(targdir, "scoscratch.ff", sep = ""))
    }
    list(scores = topKscores, inds = topKinds, guniv = guniv, 
        snpnames = rownames(inimgr@fflist[[1]]), call = theCall)
}
