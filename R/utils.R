#
#topKfeats = function(mgr, K, fn="inds1.ff", batchsize=200,
#   feat = c("score", "ind", "geneind"), ffind=1, ginds ) {
##
## given an eqtlTestsManager instance, presumably generated on one
## chromosome worth of SNP (thus ffind=1 by default), find the
## top K features per SNP in an nsnp x K ff matrix
## 
## to keep things simple, we have separate treatments of values (feat == "score") and
## names of features (feat == "geneind"), where the latter are assumed to be integer
## indices into a vector of probe names
##
## it would of course be very nice to compute the order once and apply it to both
## the scores and the names simultaneously.  i am not sure it will be beneficial
##
#     intests = mgr@fflist[[ffind]]
#     if (feat == "score") op = function(x)sort(x, decreasing=TRUE)[1:K]
##     else if (feat == "ind") op = function(x)order(x, decreasing=TRUE)[1:K]
#     else if (feat == "geneind") op = function(x)ginds[
#                                        order(x,decreasing=TRUE)[1:K]]
#     else stop("feat not recognized")
#     tmp = ffrowapply(
#       t(apply(intests[i1:i2,],1,op)),
#       X=intests, RETURN=TRUE, RETCOL=K,
#       BATCHSIZE=batchsize)
#     ff(tmp, filename=fn, dim=c(nrow(intests),K), overwrite=TRUE,
#       vmode="short")
#       }
#
#updateKfeats = function( sco1, sco2, ind1, ind2, batchsize=200 ) {
##
## will overwrite sco1, ind1 with the improvements available in sco2, ind2
##
#   snchunk = chunk(1, nrow(sco1), by=batchsize)
#   K=ncol(sco1)
#   for (i in 1:length(snchunk)) {
#      i1 = snchunk[[i]][1]
#      i2 = snchunk[[i]][2]
#      scos = cbind(sco1[i1:i2,], sco2[i1:i2,])
#      ginds = cbind(ind1[i1:i2,], ind2[i1:i2,])
#      rowwiseExtract = function(x,y) t(sapply(1:nrow(x), function(row) x[row,][y[row,]]))
#      chind = t(apply(scos, 1, function(x)order(x,decreasing=TRUE)[1:K]))
#      sco1[i1:i2,] = rowwiseExtract( scos, chind )
#      ind1[i1:i2,] = rowwiseExtract( ginds, chind )
#      }
#   invisible(NULL)
#}
#
#
#cisZero = function(mgr, snpRanges, geneRanges, radius) {
##
## this function addresses the problem of excluding same-chromosome cis tests
## up to a certain radius around each gene
##
## given a manager, all tests for SNP within radius of each gene are set to zero
##
#        dimnt = dimnames(mgr@fflist[[1]])
#        oks = gsub("rs", "", dimnt[[1]])
#        okg = dimnt[[2]]
#        srids = values(snpRanges)$RefSNP_id
#        if (!isTRUE(all(okg %in% names(geneRanges))))  {
#            warning("geneRanges does not include names/ranges for all colnames of mgr fflist")
#            okg = intersect( okg, names(geneRanges))
#            }
#        if (!isTRUE(all(oks %in% srids)))  {
#            warning("snpRanges does not include names/ranges for all rownames of mgr fflist")
# # match and numeric indexing much more efficient than names for large objects
#            oks = match(oks, srids, nomatch=0)
#            if (length(oks) == 0) stop("snpRanges does not include any snps in mgr")
#            }
# # can't use zero indices in subsetting GRanges...
#        ol = findOverlaps(snpRanges[oks[oks>0]], geneRanges[okg] + 
#            radius)
#        matm = matchMatrix(ol)
#        if (nrow(matm) > 0) {
#            mgr@fflist[[1]][matm] = 0
#        }
#    }
#
#transScores = function (smpack, snpchr = "chr1", rhs, K = 20, targdir = "tsco", 
#    geneApply = mclapply, chrnames = paste("chr", as.character(1:22), sep=""), 
#    geneRanges = NULL, snpRanges = NULL, radius = 2e+06) 
#{
##
## objective is a small-footprint accumulation of trans-eQTL tests
##  smpack is a lightweight smlSet package name
##  snpchr is the chromosome for which SNPs will be tested
##  rhs is the right hand side of formula for snp.rhs.tests in snpTests
##  K is the number of best features to be retained as we explore the transcriptome
##
##
#    if (length(chrnames) < 2) 
#        stop("must have length(chrnames) >= 2")
#    theCall = match.call()
#    require(GGtools)
##
## get an image of the expression+genotype data for SNP on specific chromosome snpchr
##
#    sms = getSS(smpack, snpchr)
#    guniv = featureNames(sms)   # universe of probes
#    smanno = gsub(".db", "", annotation(sms))
#    require(paste(smanno, ".db", sep = ""), character.only = TRUE)
#    clcnames = gsub("chr", "", chrnames)  # typical chrom nomenclature of bioconductor
#    pnameList = mget(clcnames, revmap(get(paste(smanno, "CHR", 
#        sep = ""))))
# # be sure to use only genes that are on arrays in sms
#    pnameList = lapply(pnameList, function(x) intersect(x, guniv))
#    names(pnameList) = chrnames  # now the universe of probes is split into chromsomes
#    todrop = which(sapply(pnameList, length)==0)
#    if (length(todrop)>0) pnameList = pnameList[-todrop]
#    genemap = lapply(pnameList, function(x) match(x, guniv))  # numerical indices for probes
#    nchr_genes = length(names(pnameList))
#    inimgr = eqtlTests(sms[probeId(pnameList[[chrnames[1]]]),   # start the sifting through transcriptome
#        ], rhs, targdir = targdir, runname = paste("tsc_", chrnames[1],  # testing on genes in chrom 1
#        sep = ""), geneApply = geneApply, saveSummaries = FALSE)
#    if (snpchr == chrnames[1]) {
#        if (is.null(geneRanges) || is.null(snpRanges)) 
#            stop("ranges must be supplied to exclude cis tests")
#        cisZero(inimgr, snpRanges, geneRanges, radius)   # if SNP are on chrom 1, exclude cis
#    }
#    topKinds = topKfeats(inimgr, K = K, fn = paste(targdir, "/",  # sort and save
#        snpchr, "_tsinds1_1.ff", sep = ""), feat = "geneind", 
#        ginds = genemap[[1]])
#    topKscores = topKfeats(inimgr, K = K, fn = paste(targdir, 
#        "/", snpchr, "_tssco1_1.ff", sep = ""), feat = "score", 
#        ginds = genemap[[1]])
#    unlink(filename(inimgr@fflist[[1]]))
#    for (j in 2:nchr_genes) {    # continue sifting through transcriptome
#        cat(j)
#        nxtmgr = eqtlTests(sms[probeId(pnameList[[chrnames[j]]]), 
#            ], rhs, targdir = targdir, runname = paste("tsctmp", 
#            j, sep = ""), geneApply = geneApply, saveSummaries = FALSE)
#        if (snpchr == chrnames[j]) {
#            if (is.null(geneRanges) || is.null(snpRanges)) 
#                stop("ranges must be supplied to exclude cis tests")
#            cisZero(nxtmgr, snpRanges, geneRanges, radius)
#            }
#        nxtKinds = topKfeats(nxtmgr, K = K, fn = paste(targdir, 
#            "indscratch.ff", sep = ""), feat = "geneind", ginds = genemap[[j]])
#        nxtKscores = topKfeats(nxtmgr, K = K, fn = paste(targdir, 
#            "scoscratch.ff", sep = ""), feat = "score", ginds = genemap[[j]])
#        updateKfeats(topKscores, nxtKscores, topKinds, nxtKinds)  
#        unlink(filename(nxtmgr@fflist[[1]]))   # kill off scratch materials
#        unlink(paste(targdir, "indscratch.ff", sep = ""))
#        unlink(paste(targdir, "scoscratch.ff", sep = ""))
#    }
#    list(scores = topKscores, inds = topKinds, guniv = guniv, 
#        snpnames = rownames(inimgr@fflist[[1]]), call = theCall)
#}
