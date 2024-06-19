bed_to_granges_dynamic <- function(bed,header=FALSE){
    require(GenomicRanges)
    bed <- read.delim(bed, stringsAsFactors=F, sep="\t",header=header)
    gr <- GRanges(seqnames=Rle(bed[,1]), ranges=IRanges(start=as.numeric(bed[,2]), end=as.numeric(bed[,3])))
    gr@elementMetadata@listData <- as.list(bed[4:ncol(bed)])

    return(gr)}

intersect_with_metadata <- function(gr1,gr2, ignore.strand=FALSE){
    require(GenomicRanges)
    if(ignore.strand == FALSE){ out <- gr1[findOverlaps(gr1,gr2)@from,]} else if( ignore.strand ==TRUE){ out <- gr1[findOverlaps(gr1,gr2,ignore.strand=TRUE)@from,]}
    ; return(out)}

setdiff_with_metadata <- function(gr1,gr2){ require(GenomicRanges)
                                            int <- findOverlaps(gr1,gr2); out <- gr1[setdiff(1:length(gr1),unique(int@from)),]; return(out)}

sort_gr <- function(gr){
    require(GenomicRanges)
    gr <- sort(sortSeqlevels(gr))

    return(gr)}
