setGeneric("gps2xgmml", function(x, fname='out.xgmml', names.arg, fontSize=4, col=col2hex('steelblue'), cex) standardGeneric("gps2xgmml"))
setMethod("gps2xgmml", signature(x='mds'), function(x, fname="out.xgmml", names.arg, fontSize=4, col=col2hex('steelblue'), cex) {
  if (ncol(x@points)==2) {
    xgmml2d(100*x@points, fname=fname, names.arg=names.arg, fontSize=fontSize, col=col, cex=cex)
  } else if (ncol(x@points==3)) {
    xgmml3d(100*x@points, fname=fname, names.arg=names.arg, fontSize=fontSize, col=col, cex=cex)
  }
}
)


#Export 2-dimensional MDS to XGMML
xgmml2d <- function(x, fname, names.arg, fontSize, col, cex) {
  if (missing(names.arg)) names.arg <- rownames(x)
  if (length(names.arg) != nrow(x)) stop('length(names.arg) does not match number of points to represent')
  if (length(col) == 1) col <- rep(col,nrow(x))
  if (length(col) != nrow(x)) stop('length(col) must be either 1 or match the number of points to represent')
  if (missing(cex)) cex <- 12
  if (is.numeric(cex)) cex <- as.character(cex) else stop('cex should be numeric')
  if (length(cex) == 1) cex <- rep(cex,nrow(x))
  if (length(cex) != nrow(x)) stop('length(cex) must be either 1 or match the number of points to represent')
  #if (missing(bindingSites)) bindingSites <- rep(NA,nrow(x))
  wid <- "0.1"; outline <- "#666666"
  font <- paste("SansSerif.bold-0-",fontSize,sep="")
  #Header and general description
  attrs <- list("chroGPS", "http://purl.org/dc/elements/1.1", "http://www.w3.org/1999/xlink", "http://www.w3.org/1999/02/22-rdf-syntax-ns#", "http://www.cytoscape.org", "http://www.cs.rpi.edu/XGMML")
  names(attrs) <- c('label','xmlns:dc','xmlns:xlink','xmlns:rdf','xmlns:cy','xmlns')
  out <- xmlOutputDOM(tag="graph", attrs=attrs)
  #out <- xmlOutputDOM(tag="graph", attrs=attrs, xmlDeclaration=TRUE)
  out$addTag(tag="att", attrs=list(name="documentVersion",value="1.1"))
  out$addTag(tag="att", attrs=list(name="networkMetadata"), close=FALSE)
    out$addTag(tag="rdf:RDF", close=FALSE)
      attrs <- list("http://www.cytoscape.org"); names(attrs) <- 'rdf:about'
      out$addTag(tag="rdf:Description", attrs=attrs, close=FALSE)
        out$addTag("dc", attrs=list(type="chroGPS",description="",identifier="",format="Cytoscape-XGMML"))
      out$closeTag("") #rdf:Description
    out$closeTag("") #rdf:RDF
  out$closeTag("") #networkMetaData
  #General viewing options
  out$addTag("att", attrs=list(type="string",name="backgroundColor",value="#ccccff"))
  out$addTag("att", attrs=list(type="real",name="GRAPH_VIEW_ZOOM",value="2.0"))
  out$addTag("att", attrs=list(type="real",name="GRAPH_VIEW_CENTER_X",value="0.0"))
  out$addTag("att", attrs=list(type="real",name="GRAPH_VIEW_CENTER_Y",value="0.0"))
  out$addTag("att", attrs=list(type="boolean",name="NODE_SIZE_LOCKED",value="true"))
  attrs <- list(type="string",name="__layoutAlgorithm",value="attribute-circle","true"); names(attrs)[4] <- 'cy:hidden'
  out$addTag("att", attrs=attrs)
  #Save Nodes
  nattrs <- c('type','h','w','x','y','fill','width','outline','cy:nodeTransparency','cy:nodeLabelFont','cy:nodeLabel','cy:borderLineType')
  for (i in 1:nrow(x)) {
    out$addTag("node", attrs=list(label=names.arg[i], id=i), close=FALSE)
      out$addTag("att", attrs=list(type="string",name="name",value=names.arg[i]))
      #out$addTag("att", attrs=list(type="integer",name="binding sites",value=bindingSites[i]))
      attrs <- list("ELLIPSE",cex[i],cex[i],x[i,1],x[i,2],col[i],wid,outline,"1.0",font,names.arg[i],"solid")
      names(attrs) <- nattrs
      out$addTag("graphics", attrs=attrs)
    out$closeTag("") #node
  }
  out$closeTag("graph") #graph
  saveXML(out$value(), file=fname)
}




#Export 3-dimensional MDS to XGMML
xgmml3d <- function(x, fname, names.arg, fontSize, col, cex) {
  if (missing(names.arg)) names.arg <- rownames(x)
  if (length(names.arg) != nrow(x)) stop('length(names.arg) does not match number of points to represent')
  if (length(col) == 1) col <- rep(col,nrow(x))
  if (length(col) != nrow(x)) stop('length(col) must be either 1 or match the number of points to represent')
  if (missing(cex)) cex <- 100
  if (is.numeric(cex)) cex <- as.character(cex) else stop('cex should be numeric')
  if (length(cex) == 1) cex <- rep(cex,nrow(x))
  if (length(cex) != nrow(x)) stop('length(cex) must be either 1 or match the number of points to represent')
  #if (missing(bindingSites)) bindingSites <- rep(NA,nrow(x))
  wid <- "0.1"; outline <- "#666666"
  font <- paste("SansSerif.bold-0-",fontSize,sep="")
  #Header and general description
  attrs <- list("chroGPS", "http://purl.org/dc/elements/1.1", "http://www.w3.org/1999/xlink", "http://www.w3.org/1999/02/22-rdf-syntax-ns#", "http://www.cytoscape.org", "http://www.cs.rpi.edu/XGMML")
  names(attrs) <- c('label','xmlns:dc','xmlns:xlink','xmlns:rdf','xmlns:cy','xmlns')
  out <- xmlOutputDOM(tag="graph", attrs=attrs)
  out$addTag(tag="att", attrs=list(name="documentVersion",value="1.1"))
  out$addTag(tag="att", attrs=list(name="networkMetadata"), close=FALSE)
    out$addTag(tag="rdf:RDF", close=FALSE)
      attrs <- list("http://www.cytoscape.org"); names(attrs) <- 'rdf:about'
      out$addTag(tag="rdf:Description", attrs=attrs, close=FALSE)
        out$addTag("dc", attrs=list(type="chroGPS",description="",identifier="",format="Cytoscape-XGMML"))
      out$closeTag("") #rdf:Description
    out$closeTag("") #rdf:RDF
  out$closeTag("") #networkMetaData
  #General viewing options
  out$addTag("att", attrs=list(name="name",value="chroGPS 3D",type="string"))
  out$addTag("att", attrs=list(name="shared name",value="chroGPS 3D",type="string"))
  out$addTag("graphics", close=FALSE)  
    out$addTag("att", attrs=list(name="LIGHT_AMBIENT_ALPHA",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="LIGHT_AMBIENT_COLOR",value="#ffffff",type="string"))
    out$addTag("att", attrs=list(name="LIGHT_DIFFUSE_ALPHA",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="LIGHT_DIFFUSE_COLOR",value="#ffffff",type="string"))
    out$addTag("att", attrs=list(name="LIGHT_ENABLED",value="true",type="string"))
    out$addTag("att", attrs=list(name="LIGHT_SPECULAR_COLOR",value="#ffffff",type="string"))
    out$addTag("att", attrs=list(name="LIGHT_SPECULAR_ALPHA",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="LIGHT_X_LOCATION",value="200.0",type="string"))
    out$addTag("att", attrs=list(name="LIGHT_Y_LOCATION",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="LIGHT_Z_LOCATION",value="400.0",type="string"))
    out$addTag("att", attrs=list(name="CAMERA_PITCH_ANGLE",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="CAMERA_ROLL_ANGLE",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="CAMERA_YAW_ANGLE",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="CAMERA_X_DIRECTION",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="CAMERA_Y_DIRECTION",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="CAMERA_Z_DIRECTION",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="CAMERA_X_UP",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="CAMERA_Y_UP",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="CAMERA_Z_UP",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="CAMERA_X_LOCATION",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="CAMERA_Y_LOCATION",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="CAMERA_Z_LOCATION",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="NETWORK_BACKGROUND_PAINT",value="#ccccff",type="string"))
    out$addTag("att", attrs=list(name="NETWORK_CENTER_X_LOCATION",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="NETWORK_CENTER_Y_LOCATION",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="NETWORK_CENTER_Z_LOCATION",value="0.0",type="string"))
    out$addTag("att", attrs=list(name="NETWORK_DEPTH",value="200.0",type="string"))
    out$addTag("att", attrs=list(name="NETWORK_HEIGHT",value="200.0",type="string"))
    out$addTag("att", attrs=list(name="NETWORK_SCALE_FACTOR",value="1.0",type="string"))
    out$addTag("att", attrs=list(name="NETWORK_SIZE",value="550.0",type="string"))
    out$addTag("att", attrs=list(name="NETWORK_TITLE",value="chroGPS",type="string"))
    out$addTag("att", attrs=list(name="NODE_SIZE_LOCKED",value="true",type="string"))
    out$addTag("att", attrs=list(name="SHOW_EDGE_LABELS",value="true",type="string"))
    out$addTag("att", attrs=list(name="SHOW_NODE_LABELS",value="true",type="string"))
  out$closeTag("") #graphics
  #Save Nodes
  nattrs <- c('type','x','y','z','fill','cy:nodeLabelFont','cy:nodeLabel')
  for (i in 1:nrow(x)) {
    out$addTag("node", attrs=list(label=names.arg[i], id=i), close=FALSE)
      out$addTag("att", attrs=list(type="string",name="name",value=names.arg[i]))
      attrs <- list("ELLIPSE",x[i,1],x[i,2],x[i,3],col[i],font,names.arg[i])
      names(attrs) <- nattrs
      out$addTag("graphics", attrs=attrs, close=FALSE)
        out$addTag("att", attrs=list(name="NODE_HEIGHT",value=cex[i],type="string"))
        out$addTag("att", attrs=list(name="NODE_WIDTH",value=cex[i],type="string"))
        out$addTag("att", attrs=list(name="NODE_DEPTH",value=cex[i],type="string"))
      out$closeTag("") #graphics
    out$closeTag("") #node
  }
  out$closeTag("graph") #graph
  saveXML(out$value(), file=fname)
}
