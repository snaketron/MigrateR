
get_pairs <- function(x) {
    p <- extract(x$f, par = "mu_group")$mu_group
    if(ncol(p)==1) {
        stop("only one treatment group: nothing to compare")
    }
    
    ds <- vector(mode = "list", length = ncol(p)*ncol(p))
    ct <- 1
    for(i in 1:ncol(p)) {
        for(j in 1:ncol(p)) {
            d <- p[,i]-p[,j]
            d_M <- mean(d)
            d_HDI <- get_hdi(vec = d, hdi_level = 0.95)
            pmax <- get_pmax(d)
            ds[[ct]] <- data.frame(group_id_x = i, 
                                   group_id_y = j, 
                                   delta_M = d_M,
                                   delta_L95 = d_HDI[1],
                                   delta_H95 = d_HDI[2],
                                   pmax = pmax)
            ct <- ct + 1
        }
    }
    ds <- do.call(rbind, ds)
    
    ds <- merge(x = ds, y = x$s$mu_group[, c("group_id", "group")],
                    by.x = "group_id_x", by.y = "group_id")
    ds <- merge(x = ds, y = x$s$mu_group[, c("group_id", "group")],
                    by.x = "group_id_y", by.y = "group_id")
    colnames(ds) <- gsub(pattern = "\\.", replacement = '_', x = colnames(ds))
    
    g <- ggplot(data = ds)+
        geom_tile(aes(x = group_x, y = group_y), col = "white", fill = "#eeeeee")+
        geom_point(aes(x = group_x, y = group_y, size = pmax, col = delta_M))+
        geom_text(aes(x = group_x, y = group_y, 
                      label = round(x = pmax, digits = 2)), size = 2)+
        scale_color_distiller(name = expression(delta), palette = "Spectral")+
        scale_radius(name = expression(pi), limits = c(0, 1))+
        theme_bw(base_size = 10)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        xlab(label = '')+
        ylab(label = '')
    
    g
    
    return(list(ds = ds, plot = g))
}


get_pmax <- function(x) {
  if(all(x==0)) {
    return(0)
  }
  l <- length(x)
  return(2*max(sum(x<0)/l, sum(x>0)/l)-1)
}

# Description:
# Computes HDI for vector vec and hdi_level (e.g. 0.95)
# Taken (and renamed) from "Doing Bayesian Analysis", section 25.2.3 R code
# for computing HDI of a MCMC sample
get_hdi <- function(vec, hdi_level) {
  sortedPts <- sort(vec)
  ciIdxInc <- floor(hdi_level * length(sortedPts))
  nCIs = length(sortedPts) - ciIdxInc
  ciWidth = rep(0 , nCIs)
  for (i in 1:nCIs) {
    ciWidth[i] = sortedPts[i + ciIdxInc] - sortedPts[i]
  }
  HDImin = sortedPts[which.min(ciWidth)]
  HDImax = sortedPts[which.min(ciWidth) + ciIdxInc]
  HDIlim = c(HDImin, HDImax)
  return(HDIlim)
}

