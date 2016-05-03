# Routines to retrieve data from modEncode, both in a generic way and directly using their API (below)

getURL <- function(urls,filenames,extension='.gff3',method='internal')
  {
    if (length(urls)!=length(filenames)) stop('Filenames not matching urls length')
    for (i in 1:length(urls))
        {
          out <- paste(filenames[i],extension,sep='')
          print(sprintf('*** Downloading %s from %s ***',out,urls[i]))
          download.file(urls[i],out,method=method,quiet=TRUE)
        }
    print(sprintf('*** Downloaded specified URLs in %s',getwd()))    
  }

# Deprecated function, too complex and too subject to changes from the modEncode query system
getmodEncodeBS <- function(type='query',ids,release='21',organism="Drosophila+melanogaster",celltype="S2-DRSC",exptype="ChIP-chip",track='gff3',method="internal")
{
  if (type=='query')
    {
      # First we download the table of results from which to build the individual dCCids to download
      # Base URL
      urlroot <- paste("http://intermine.modencode.org/release-",release,"/service/query/results?query=",sep='')
      # Building query
      query <- paste('%3Cquery+name%3D%22%22+model%3D%22genomic%22+view%3D%22BindingSite.submissions.title+BindingSite.submissions.DCCid+BindingSite.submissions.design+BindingSite.submissions.experimentDate+BindingSite.submissions.embargoDate%22+sortOrder%3D%22BindingSite.submissions.title+asc%22+constraintLogic%3D%22A+and+C+and+D%22%3E%3Cconstraint+path%3D%22BindingSite.organism%22+code%3D%22A%22+op%3D%22LOOKUP%22+value%3D%22',organism,'%22+extraValue%3D%22%22%2F%3E%3Cconstraint+path%3D%22BindingSite.submissions.experimentType%22+code%3D%22C%22+op%3D%22LIKE%22+value%3D%22',exptype,'%22%2F%3E%3Cconstraint+path%3D%22BindingSite.submissions.properties%22+type%3D%22CellLine%22%2F%3E%3Cconstraint+path%3D%22BindingSite.submissions.properties.name%22+code%3D%22D%22+op%3D%22%3D%22+value%3D%22',celltype,'%22%2F%3E%3C%2Fquery%3E',sep='')
      # URL tail, format for table
      urltail <- "&format=csv"
      # Full URL
      fullurl=paste(urlroot,query,urltail,sep='')
      format <- "csv" # Format for table of contents to download
      print(sprintf('*** Downloading modEncode release %s GFF Binding Sites for %s, %s, %s ***',release,organism,celltype,exptype))
      download.file(fullurl,paste('result',format,sep='.'),method=method,quiet=TRUE) # Method will depend on system?
      # Now proceed with all downloads
      # Define base URL (will depend on data release)
      urlroot <- paste("http://intermine.modencode.org/release-",release,"/",sep='')
      urlbody <- paste("features.do?type=submission&action=export&format=",track,"&submission=",sep='')
      urltail="&feature=BindingSite"
      # Read table of contents
      toc <- read.csv(paste('result',format,sep='.'),header=FALSE,as.is=TRUE); colnames(toc) <- c('Name','ID','Type','ReleaseDate','EmbargoDate')
      toc$Name <- gsub(' ','-',toc$Name)
      for (i in 1:nrow(toc))
        {
          fullurl <- paste(urlroot,urlbody,toc[i,'ID'],urltail,sep='')
          #out <- paste("'",toc[i,'ID'],'-',toc[i,'Name'],'.',track,"'",sep='')
          out <- paste(toc[i,'ID'],'-',toc[i,'Name'],'.',track,sep='')
          print(sprintf('*** Downloading BindingSite information for %s:%s ***',toc[i,'ID'],toc[i,'Name']))
          download.file(fullurl,out,method=method,quiet=TRUE)
          #url.show(fullurl,file=out,method=method)
        }
      print(sprintf('*** Downloaded modEncode BindingSite information in %s',getwd()))
      return(toc)
    }
  else if (type=='batch')
    # We have a list of modencode identifiers
    {
      if (length(ids)<1) stop('Invalid modEncode identifier list')
      #ids <- paste('modENCODE',ids,sep='_')
      # Now proceed with all downloads
      # Define base URL (will depend on data release)
      urlroot <- paste("http://intermine.modencode.org/release-",release,"/",sep='')
      urlbody <- paste("features.do?type=submission&action=export&format=",track,"&submission=",sep='')
      urltail="&feature=BindingSite"
      for (id in ids)
        {
          fullurl <- paste(urlroot,urlbody,id,urltail,sep='')
          out <- paste(id,'.',track,sep='')
          print(sprintf('*** Downloading BindingSite information for %s ***',id))
          download.file(fullurl,out,method=method,quiet=TRUE)
          #url.show(fullurl,file=out,method=method)
        }
      print(sprintf('*** Downloaded modEncode BindingSite information in %s',getwd()))
      return(ids)
    }
}
