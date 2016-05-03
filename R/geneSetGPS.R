# Given a list of FlyBase IDs, returns the coordinates of those genes in a chroGPS-genes map

# Methods & Functions
# X is the FBID x variables table used to generate the map
setGeneric("geneSetGPS", function(x, m, genes, uniqueCount=TRUE,...) standardGeneric("geneSetGPS"))

setMethod("geneSetGPS", signature(x="data.frame",m="mds",genes="character"), function(x, m, genes, uniqueCount=TRUE, ...) {
  if (uniqueCount) # We have a mapping of FBIDs to unique genes according to their variables
    {
      #xu <- uniqueCount(x)
      #rownames(xu) <- xu$u
      txt <- paste("x$u <- paste(",paste("x[,",1:ncol(x),"]",collapse=","),", sep=',')",sep="")
      eval(parse(text=txt))
      genes <- genes[genes %in% rownames(x)]
      genesPattern <- x[genes,'u']; names(genesPattern) <- genes
      #xu <- xu[rownames(m@points),] # Reindex to revert splitMDS reshuffling if necessary
      #sel <- xu$u %in% genesPattern
      genesPattern <- genesPattern[genesPattern %in% rownames(m@points)]
      selpoints <- m@points[genesPattern,]
      rownames(selpoints) <- names(genesPattern)
      new("mds",points=selpoints,R.square=0)
    }
  else # We can use the original mappings
    {
      new("mds",points=m@points[genes,],R.square=0)
    }
}
)

setMethod("geneSetGPS", signature(x="matrix",m="mds",genes="character"), function(x, m, genes, uniqueCount=TRUE,...) {
  geneSetGPS(x=as.data.frame(x), m, genes, uniqueCount,...)
}
)
