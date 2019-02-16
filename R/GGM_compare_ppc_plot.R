#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
GGM_compare_ppc_plot <- function(X, alpha, bw, log_transform = NULL){

  # change object name
  fit <- X

  if(class(fit) != "GGM_compare_pcc") stop("incorrect object class")

  # number of contrasts
  k <- length(fit)

  # predicitve risk for k contrasts
  risk <- unlist(lapply(1:k, function(x) t(fit[[x]]$jsd_scores)))

  # contrast names
  cont_names <- unlist(lapply(1:k, function(x) fit[[x]]$contrast))
  temp <- stringr::str_remove_all(cont_names, "_")

  # reformat names
  cont_names <- BGGM:::axis_ticks_helper(temp)

  # observed JSD
  obs_jsd <- unlist(lapply(1:k, function(x) fit[[x]]$obs_jsd) )

  # data frame for plotting the observed JSD
  dat_point <- cbind.data.frame(contrast_names = unique(cont_names),
                                obs_jsd = obs_jsd)


  # contrast names for plotting
  contrast_names <- rep(cont_names, each = length(risk) / k)

  # observed JSD for plotting
  obs_jsd <- rep(obs_jsd, each = length(risk) / k)

  # combine for plotting
  dat_plt <- cbind.data.frame(risk, contrast_names, obs_jsd)

  # not transformed
  if(is.null(log_transform)){

    plt <- ggplot(dat_plt, aes(group = contrast_names,
                               risk,
                               y = contrast_names)) +
      # ridges
      geom_density_ridges(aes(fill = contrast_names),
                          alpha = alpha,
                          show.legend = F,
                          bandwidth = bw) +
      # observed JSD
      geom_point(data = dat_point, aes(x = obs_jsd,
                                       group = contrast_names)) +
      ylab("Constrasts")
  }

  if(isTRUE(log_transform)){

    plt <- ggplot(dat_plt, aes(group = contrast_names,
                               log(risk),
                               y = contrast_names)) +
      # ridges
      geom_density_ridges(aes(fill = contrast_names),
                          alpha = alpha,
                          show.legend = F,
                          bandwidth = bw) +
      # observed JSD
      geom_point(data = dat_point, aes(x = log(obs_jsd),
                                       group = contrast_names)) +
      ylab("Constrasts")
    }
  plt
}

