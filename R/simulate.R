
sim <- function(N_well_cells,
                N_plate, 
                N_group, 
                N_well_reps, 
                shape, 
                sigma_bplate, 
                sigma_wplate, 
                alpha_plate, 
                mu_group) {
    
    get_meta <- function(data_in) {
        meta <- c()
        well_index <- 1
        for(plate_index in 1:data_in$N_plate) {
            for(group_index in 1:data_in$N_group) {
                for(w in 1:data_in$N_well_reps) {
                    r <- data.frame(sample = paste0("s", well_index),
                                    well_id = well_index, 
                                    group = paste0("g", group_index),
                                    group_id = group_index,
                                    plate = paste0("p", plate_index),
                                    plate_id = plate_index)
                    meta <- rbind(meta, r)
                    well_index = well_index + 1;
                }
            }
        }
        return(meta)
    }
    
    data_in <- list(N_plate = N_plate,
                    N_group = N_group,
                    N_well_reps = N_well_reps,
                    N_well_cells = N_well_cells,
                    shape = shape,
                    sigma_bplate = sigma_bplate,
                    sigma_wplate = sigma_wplate,
                    alpha_plate = alpha_plate,
                    mu_group = mu_group)
    
    message("simulation... \n")
    
    # sim data from model
    f <- sampling(object = stanmodels$S,
                  data = data_in,
                  chains = 1, 
                  cores = 1,
                  warmup = 1,
                  iter = N_well_cells+1,
                  algorithm = "Fixed_param")
    
    # prep data
    y <- extract(object = f, par = "y_hat_sample")$y_hat_sample
    y <- melt(y)
    colnames(y) <- c("iteration", "well_id", "v")
    y$iteration <- NULL
    
    # meta
    m <- get_meta(data_in = data_in)
    y <- merge(x = y, y = m, by.x = "well_id", by.y = "well_id", all.x = TRUE)
    
    return(y)
}
