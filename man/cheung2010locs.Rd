\name{cheung2010locs}
\alias{cheung2010locs}
\docType{data}
\title{
cheung2010locs gives locations for HapMap phase 3 SNPs for CEPH CEU
}
\description{
cheung2010locs gives locations for HapMap phase 3 SNPs for CEPH CEU
}
\usage{data(cheung2010locs)}
\format{
  The format is:\cr
Formal class 'GRangesList' [package "GenomicRanges"] with 5 slots\cr
  ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots\cr
  .. .. ..@ end            : int [1:22] 112259 226087 320873 406401 493642 585295 660727 735268 798907 871895 ...\cr
  .. .. ..@ NAMES          : chr [1:22] "chr1" "chr2" "chr3" "chr4" ...\cr
  .. .. ..@ elementType    : chr "integer"\cr
  .. .. ..@ elementMetadata: NULL\cr
  .. .. ..@ metadata       : list()\cr
  ..@ unlistData     :Formal class 'GRanges' [package "GenomicRanges"] with 6 slots\cr
  .. .. ..@ seqnames       :Formal class 'Rle' [package "IRanges"] with 4 slots\cr
  .. .. .. .. ..@ values         : Factor w/ 22 levels "chr1","chr2",..: 1 2 3 4 5 6 7 8 9 10 ...\cr
  .. .. .. .. ..@ lengths        : int [1:22] 112259 113828 94786 85528 87241 91653 75432 74541 63639 72988 ...\cr
  .. .. .. .. ..@ elementMetadata: NULL\cr
  .. .. .. .. ..@ metadata       : list()\cr
  .. .. ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots\cr
  .. .. .. .. ..@ start          : int [1:1365825] 554340 557616 711153 742429 742584 743268 744045 750775 751010 751595 ...\cr
  .. .. .. .. ..@ width          : int [1:1365825] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. .. .. .. ..@ NAMES          : chr [1:1365825] "rs6650104" "rs11510103" "rs12565286" "rs3094315" ...\cr
  .. .. .. .. ..@ elementType    : chr "integer"\cr
  .. .. .. .. ..@ elementMetadata: NULL\cr
  .. .. .. .. ..@ metadata       : list()\cr
  .. .. ..@ strand         :Formal class 'Rle' [package "IRanges"] with 4 slots\cr
  .. .. .. .. ..@ values         : Factor w/ 3 levels "+","-","*": 1\cr
  .. .. .. .. ..@ lengths        : int 1365825\cr
  .. .. .. .. ..@ elementMetadata: NULL\cr
  .. .. .. .. ..@ metadata       : list()\cr
  .. .. ..@ elementMetadata:Formal class 'DataFrame' [package "IRanges"] with 6 slots\cr
  .. .. .. .. ..@ rownames       : NULL\cr
  .. .. .. .. ..@ nrows          : int 1365825\cr
  .. .. .. .. ..@ listData       :List of 1\cr
  .. .. .. .. .. ..$ assignment: Factor w/ 7 levels "./.","A/C","A/G",..: 3 3 5 3 3 2 3 6 6 6 ...\cr
  .. .. .. .. ..@ elementType    : chr "ANY"\cr
  .. .. .. .. ..@ elementMetadata: NULL\cr
  .. .. .. .. ..@ metadata       : list()\cr
  .. .. ..@ seqinfo        :Formal class 'Seqinfo' [package "GenomicRanges"] with 3 slots\cr
  .. .. .. .. ..@ seqnames   : chr [1:22] "chr1" "chr2" "chr3" "chr4" ...\cr
  .. .. .. .. ..@ seqlengths : int [1:22] NA NA NA NA NA NA NA NA NA NA ...\cr
  .. .. .. .. ..@ is_circular: logi [1:22] NA NA NA NA NA NA ...\cr
  .. .. ..@ metadata       : list()\cr
  ..@ elementType    : chr "GRanges"\cr
  ..@ elementMetadata:Formal class 'DataFrame' [package "IRanges"] with 6 slots\cr
  .. .. ..@ rownames       : NULL\cr
  .. .. ..@ nrows          : int 22\cr
  .. .. ..@ listData       : Named list()\cr
  .. .. ..@ elementType    : chr "ANY"\cr
  .. .. ..@ elementMetadata: NULL\cr
  .. .. ..@ metadata       : list()\cr
  ..@ metadata       : list()\cr
}
\details{
obtained as the 'support' elements after \code{\link[chopsticks]{read.HapMap.data}} was run
on the files noted below.
}
\source{
for example, 

\url{http://hapmap.ncbi.nlm.nih.gov/downloads/genotypes/hapmap3_r3/hapmap_format/polymorphic/genotypes_chr17_CEU_phase3.3_nr.b36_fwd.txt.gz}

}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(cheung2010locs)
names(cheung2010locs)
cheung2010locs[[1]][1:3]
}
\keyword{datasets}
