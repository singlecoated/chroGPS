# package caTools doesn't export the combs function correctly in all its versions, so we copy it directly here
# see info from combs function from caTools package for more info.
# combs function from caTools version 1.12
combs <- function (v, k) 
{
    n = length(v)
    if (n == k) 
        P = matrix(v, 1, n)
    else if (k == 1) 
        P = matrix(v, n, 1)
    else if (k == n - 1) 
        P = matrix(rep(v, each = n - 1), n, n - 1)
    else if (k < n) {
        P = matrix(0, 0, k)
        if (k < n & k > 1) {
            for (i in 1:(n - k + 1)) {
                Q = combs(v[(i + 1):n], k - 1)
                j = nrow(Q)
                P = rbind(P, cbind(rep(v[i], j), Q))
            }
        }
    }
    else stop("combs: number m has to be smaller or equal to length of vector v")
    return(P)
}

domainDist <- function(d,gps='factors',domain,type='intra',col='white',avg=FALSE,plot=TRUE, ...)
# Returns intra or interdomain distances, d is a distance matrix as returned by distGPS.
  {
    if (gps=='factors')
      # Domain is a character vector with factor (protein) names
      {
        if (length(domain)!=mean(nrow(d),ncol(d))) stop('Distance elements not matching domain length or not square distance matrix')
        if (type=='inter')
          {
            domains <- combs(unique(domain),k=2)
            dist.domain <- vector('list',nrow(domains))
            dist.domain <- apply(domains, 1, function(x,d,domain)
                                 {
                                   sel1 <- x[1]==domain
                                   sel2 <- x[2]==domain
                                   res <- d[sel1,sel2]
                                 },d,domain )
            names(dist.domain) <- apply(domains,1,paste,collapse='.')
            labs <- names(dist.domain)
            if (avg) { dist.domain$Avg <- unlist(dist.domain[1:length(dist.domain)]); col <- c(col,'white'); labs <- c(labs,'Average') }
            if (plot) boxplot(dist.domain[1:length(dist.domain)],names=labs,xlab='Domain',ylab='Distance',ylim=c(0,1),...)
          }
        else if (type=='intra')
          {
            labs <- sort(unique(domain))
            col <- unique(col)
            dist.domain <- lapply(sort(unique(domain)), function(x,d,domain)
                                  {
                                    sel <- x==domain
                                    res <- d[sel,sel]
                                    res <- res[upper.tri(res)]
                                  },d,domain )
            names(dist.domain) <- sort(unique(domain))
            if (avg) { dist.domain$Avg <- unlist(dist.domain[1:length(dist.domain)]); col <- c(col,'white'); labs <- c(labs,'Average') }
            if (plot) boxplot(dist.domain[1:length(dist.domain)],names=labs,xlab='Domain',ylab='Distance',col=col,ylim=c(0,1),...)
          }
        return(dist.domain)
      }
    else if (gps=='genes')
      # Domain is a named list (names=factors) defining distance matrix positions (genes) having that factor mark
      {
        if (type=='inter')
          {
            labs <- names(domain)
            col <- unique(col)
            domains <- combs(names(domain),k=2)
            dist.domain <- vector('list',nrow(domains))
            dist.domain <- apply(domains, 1, function(x,d,domain)
                                 {
                                   res <- d[domain[[x[1]]],domain[[x[2]]]]
                                   res <- res[upper.tri(res)]
                                   return(res)
                                 },d,domain )
            names(dist.domain) <- apply(domains,1,paste,collapse='\n')
            labs <- names(dist.domain)
            if (avg) { dist.domain$Avg <- unlist(dist.domain[1:length(dist.domain)]); col<- c(col,'white'); labs <- c(labs,'Average') }
            if (plot) boxplot(dist.domain[1:length(dist.domain)],ylim=range(dist.domain[1:length(dist.domain)]),names=labs,xlab='Domain',ylab='Distance',col=col,ylim=c(0,1),outline=FALSE,...)
          }
        else if (type=='intra')
          {
            labs <- names(domain)
            col <- unique(col)
            dist.domain <- lapply(domain, function(domain,d)
                                  {
                                    res <- d[domain,domain]
                                    res <- res[upper.tri(res)]
                                    return(res)
                                  },d)
            labs <- names(dist.domain)
            if (avg) { dist.domain$Avg <- unlist(dist.domain[1:length(dist.domain)]); col <- c(col,'white'); labs <- c(labs,'Average') }
            if (plot) boxplot(dist.domain[1:length(dist.domain)],ylim=range(dist.domain[1:length(dist.domain)]),names=labs,xlab='Domain',ylab='Distance',col=col,ylim=c(0,1),outline=FALSE,...)
          }
        return(dist.domain)
      }
  }
