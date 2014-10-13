# Functions to retrieve data from GFF3 files (modEncode) and storing them into a RangedDataList object

getAttributeField <-
function (x, field, attrsep = ";") {
     s = strsplit(x, split = attrsep, fixed = TRUE)
     sapply(s, function(atts) {
         a = strsplit(atts, split = "=", fixed = TRUE)
         m = match(field, sapply(a, "[", 1))
         if (!is.na(m)) {
             rv = a[[m]][2]
         }
         else {
             rv = as.character(NA)
         }
         return(rv)
     })
}

gffRead <-
function(gffFile, nrows = -1) {
     cat("Reading ", gffFile, ": ", sep="")
     gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
     header=FALSE, comment.char="#", nrows = nrows,
     colClasses=c("character", "character", "character", "integer",  
"integer",
     "character", "character", "character", "character"))
     colnames(gff) = c("seqname", "source", "feature", "start", "end",
             "score", "strand", "frame", "attributes")
     cat("found", nrow(gff), "rows with classes:",
         paste(sapply(gff, class), collapse=", "), "\n")
     stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
     return(gff)
}

gff2RDList <-
function(filenames,listnames=NULL,dir,quote=NULL,chrprefix='')
  {
    if (is.null(listnames)) listnames <- filenames
    res.pos <- res.neg <- vector('list',length(filenames))
    names(res.pos) <- names(res.neg) <- listnames
    gfffiles <- filenames
    if (!is.null(quote)) gfffiles <- paste(quote,gfffiles,quote,sep='')
    for (i in 1:length(res.pos))
      {
        gff <- gffRead(file.path(dir,gfffiles[i])) # Read GFF file
        gff.pos <-  subset(gff,score>0) # To keep enriched regions only
        gff.neg <-  subset(gff,score<0) # To keep depleted regions only
        rd.pos <- RangedData(IRanges(start=gff.pos$start,end=gff.pos$end),space=paste(chrprefix,gff.pos$seqname,sep=''),score=as.numeric(gff.pos$score),ID=getAttributeField(gff.pos$attributes,'ID'))
        rd.neg <- RangedData(IRanges(start=gff.neg$start,end=gff.neg$end),space=paste(chrprefix,gff.neg$seqname,sep=''),score=as.numeric(gff.neg$score),ID=getAttributeField(gff.neg$attributes,'ID'))
        res.pos[[i]] <- rd.pos
        res.neg[[i]] <- rd.neg
      }
    res.pos <- RangedDataList(res.pos)
    res.neg <- RangedDataList(res.neg)
    return(list(Enriched=res.pos,Depleted=res.neg))
  }

