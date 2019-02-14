#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
GGM_compare_ppc_plot <- function(x){

  if(class(x) != "GGM_compare_pcc")stop("incorrect object class")

  k <- length(x)
  risk <- unlist(lapply(1:k, function(x) t(fit[[x]]$jsd_scores)) )

  cont_names <- unlist(lapply(1:k, function(x) fit[[x]]$contrast) )
  obs_jsd <- unlist(lapply(1:k, function(x) fit[[x]]$obs_jsd) )

  cont_names <- gsub("  ", "", gsub("_", " ", gsub("Y", " ", cont_names, fixed=TRUE)))



  dat_point <- cbind.data.frame(contrast_names = unique(cont_names), obs_jsd)


  contrast_names <- rep(cont_names, each = length(risk) / k)
  obs_jsd <- rep(obs_jsd, each = length(risk) / k)


  dat_plt <- cbind.data.frame(risk, contrast_names, obs_jsd)

  ggplot(dat_plt, aes(x = risk)) +
    geom_density(show.legend = F, aes(fill = contrast_names),  color = NA,alpha = 0.85) +
    geom_point(data = dat_point, aes(x = obs_jsd, y = 0)) +
    facet_wrap(~ contrast_names)



  plt <- ggplot(dat_plt, aes(group = contrast_names, log(risk), y= contrast_names)) +
    geom_density_ridges(aes(fill = contrast_names), alpha = 0.75, show.legend = F,  bandwidth = .1) +

    geom_point(data = dat_point, aes(x = log(obs_jsd), group = contrast_names)) +
    ylab("Constrasts")

  plt

}
