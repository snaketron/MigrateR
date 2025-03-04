
get_ppc <- function(x) {
    e <- rstan::extract(object = x$f, par = "y_hat_sample")$y_hat_sample
    e <- reshape2::melt(data = e)
    colnames(e) <- c("iter", "well_id", "yhat")
    q <- x$x$d[, c("well_id", "compound", "dose", "group", "plate")]
    q <- q[duplicated(q)==F,]
    e <- merge(x = e, y = q, all.x = T)
    
    g <- ggplot()+
        facet_grid(compound~plate, scales = "free")+
        geom_violin(data = e, aes(x = as.factor(dose), 
                                  y = yhat, linetype = plate), 
                    fill = NA, col = "red", linetype = "dashed")+
        geom_sina(data = x$x$d, aes(x = as.factor(dose), 
                                    y = sv, col = plate), 
                  col = "black", size = 0.3)+
        theme_bw(base_size = 10)+
        theme(legend.position = "none")+
        theme(strip.text.x = element_text(
            margin = margin(0.02,0,0.02,0, "cm")))+
        xlab(label = "dose")+
        scale_y_continuous(name = "s", breaks = scales::pretty_breaks(3))
    return(g)
}
