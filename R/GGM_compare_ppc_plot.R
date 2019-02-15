#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
GGM_compare_ppc_plot <- function(X, alpha, bw, log_transform = NULL){

  fit <- X

  if(class(x) != "GGM_compare_pcc")stop("incorrect object class")

  k <- length(fit)
  risk <- unlist(lapply(1:k, function(x) t(fit[[x]]$jsd_scores)) )

  cont_names <- unlist(lapply(1:k, function(x) fit[[x]]$contrast) )
  obs_jsd <- unlist(lapply(1:k, function(x) fit[[x]]$obs_jsd) )

  temp <- stringr::str_remove_all(cont_names, "_")


  cont_names <- BGGM:::axis_ticks_helper(temp)



  dat_point <- cbind.data.frame(contrast_names = unique(cont_names), obs_jsd)


  contrast_names <- rep(cont_names, each = length(risk) / k)
  obs_jsd <- rep(obs_jsd, each = length(risk) / k)


  dat_plt <- cbind.data.frame(risk, contrast_names, obs_jsd)


  if(is.null(log_transform)){

    plt <- ggplot(dat_plt, aes(group = contrast_names,
                               risk,
                               y = contrast_names)) +
      geom_density_ridges(aes(fill = contrast_names),
                          alpha = alpha,
                          show.legend = F,
                          bandwidth = bw) +

      geom_point(data = dat_point, aes(x = obs_jsd,
                                       group = contrast_names)) +
      ylab("Constrasts")
  }

  if(isTRUE(log_transform)){

    plt <- ggplot(dat_plt, aes(group = contrast_names,
                               log(risk),
                               y = contrast_names)) +
      geom_density_ridges(aes(fill = contrast_names),
                          alpha = alpha,
                          show.legend = F,
                          bandwidth = bw) +

      geom_point(data = dat_point, aes(x = log(obs_jsd),
                                       group = contrast_names)) +
      ylab("Constrasts")



  }


  plt

}

