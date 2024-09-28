

get_boot_profiles <- function(x, gs, hc_dist, hc_link, main_ph) {
  e <- extract(x$f, par = "eff_group_mu")$eff_group_mu[, gs]
  e <- e[sample(x = 1:nrow(e), size = min(nrow(e), 1000), replace = FALSE),]

  meta <- o$s$eff_group_mu[, c("g", "treatment", "dose")]
  meta <- meta[order(meta$g, decreasing = F),]
  meta <- meta[meta$g %in% gs, ]
  meta$g <- NULL

  boot_ph <- c()
  for(i in 1:nrow(e)) {
    u <- data.frame(g = 1:ncol(e), mu = e[i, ])
    u <- cbind(u, meta)

    q <- acast(data = u, formula = treatment~dose, value.var = "mu")

    # hclust
    hc <- hclust(dist(q, method = hc_dist), method = hc_link)
    ph <- as.phylo(x = hc)

    if(i == 1) {
      boot_ph <- ph
    }
    else {
      boot_ph <- c(boot_ph, ph)
    }
  }
  clades <- prop.clades(phy = main_ph,
                        x = boot_ph,
                        part = NULL,
                        rooted = is.rooted(main_ph))

  # add bootstrap
  main_ph$node.label <- clades

  # b = 0 for these nodes
  na_nodes <- which(is.na(main_ph$node.label))
  if(length(na_nodes)!=0) {
    main_ph$node.label[na_nodes] <- 0
  }
  return(list(main_ph = main_ph, boot_ph = boot_ph))
}


get_profiles <- function(x,
                         hc_link = "average",
                         hc_dist = "euclidean",
                         select_ds,
                         select_ts) {
  eg <- x$s$eff_group_mu
  es <- x$s$eff_sample

  if(missing(select_ds)==FALSE) {
    if(any(!select_ds %in% unique(eg$dose))) {
      stop("selected doses not found in data")
    }
    eg <- eg[eg$dose %in% select_ds, ]
    es <- es[es$dose %in% select_ds, ]
  }

  if(missing(select_ts)==FALSE) {
    if(any(!select_ts %in% unique(eg$treatment))) {
      stop("selected doses not found in data")
    }
    eg <- eg[eg$treatment %in% select_ts, ]
    es <- es[es$treatment %in% select_ts, ]
  }

  q <- acast(data = eg, formula = treatment~dose, value.var = "mean")

  # hclust
  hc <- hclust(dist(q, method = hc_dist), method = hc_link)
  ph <- as.phylo(x = hc)


  bt <- get_boot_profiles(x = x, gs = eg$g, hc_dist = hc_dist,
                          hc_link = hc_link, main_ph = ph)


  tree <- ggtree(bt$main, linetype='solid')+
    geom_point2(mapping = aes(subset=isTip==FALSE),size = 0.5, col = "black")+
    geom_tippoint(size = 2, fill = "white", shape = 21)+
    geom_tiplab(color='black', as_ylab = T, align = TRUE)+
    layout_rectangular()+
    theme_bw(base_size = 10)+
    scale_x_continuous(labels = abs)+
    geom_nodelab(geom='text', color = "#4c4c4c" ,size = 2.75, hjust=-0.2,
                 mapping = aes(label=label,subset=isTip==FALSE))

  tree <- revts(tree)

  t <- tree$data
  t <- t[order(t$y, decreasing = FALSE), ]
  tips <- t$label[t$isTip==TRUE]

  q <- eg
  q$treatment <- factor(q$treatment, levels = rev(tips))

  g <- ggplot(data = q)+
    facet_grid(treatment~., switch = "y")+
    geom_hline(yintercept = 0, linetype = "dashed", col = "gray")+
    geom_point(aes(x = dose, y = mean))+
    geom_errorbar(aes(x = dose, y = mean, ymin = X2.5., ymax = X97.5.), width = 0)+
    scale_y_continuous(position = "right", breaks = scales::pretty_breaks(n = 5))+
    theme_bw(base_size = 10)+
    theme(strip.text.y = element_text(margin = margin(0.01,0.01,0.01,0.01, "cm")))


  q <- es[es$treatment %in% q$treatment, ]
  q$treatment <- factor(q$treatment, levels = rev(tips))
  g2 <- ggplot(data = q)+
    facet_wrap(facets = treatment~replicate, nrow = length(unique(q$treatment)), switch = "y", scales = "free_y")+
    geom_errorbar(aes(x = dose, y = mean, ymin = X2.5., ymax = X97.5.), width = 0, alpha = 0.5)+
    geom_line(aes(x = dose, y = mean))+
    geom_point(aes(x = dose, y = mean))+
    scale_y_continuous(position = "right", breaks = scales::pretty_breaks(n = 3))+
    theme_bw(base_size = 10)+
    theme(legend.position = "none", strip.text.y = element_text(margin = margin(0.01,0.01,0.01,0.01, "cm")))


  gs <- (tree|g|g2)+
    plot_annotation(tag_levels = 'A')

  gs
  return(gs)
}
