#' @name print.BGGM
#' @title  Print Method for \code{BGGM} Objects
#'
#' @param x An object of class \code{BGGM}.
#'
#' @param ... currently ignored
#'
#' @importFrom methods is
#' @export

print.BGGM <- function(x, ...) {
  # print estimate methods

  if(is(x, "estimate")) {

    if(is(x, "default")) {

      print_estimate(x, ...)

    } else if (is(x, "fitted")) {

      print_fitted(x, ...)

    } else if (is(x, "predict")) {
      print_predict(x, ...)

    } else if (is(x, "map")) {

      print_map(x, ...)

    } else if (is(x, "coef")) {

      print_coef(x, ...)

    } else if (is(x, "ggm_compare_estimate")) {

      if (is(x, "summary")) {

        print_summary_ggm_estimate_compare(x, ...)

      } else {

        print_ggm_compare(x, ...)

      }

      } else if (is(x, "ggm_compare_ppc")) {

        print_ggm_compare_ppc(x, ...)

    } else if (is(x, "metric")) {

      if (is(x, "summary")) {

        print_summary_metric(x)

      } else {

        if (is(x, "R2")) {

          print_summary_metric(summary(x, ...))

        } else {

          print_summary_metric(summary(x, ...))

        }
      }
      # end metric
    }  else if (is(x, "summary.estimate")) {

      print_summary_estimate(x)

    } else if (is(x, "post.pred")) {

      print_post_pred(x, ...)

    }  else if (is(x, "select.estimate")) {

      print_select_estimate(x, ...)

    } else if (is(x, "select.ggm_compare_estimate")) {

      print_select_ggm_compare_estimate(x, ...)

    }

  }
  # explore methods
  if (is(x, "explore")) {

    if(is(x, "default")){

      print_explore(x, ..)
    }

    if(is(x, "summary_explore")){

      print_summary_estimate(x,...)

    }

    if (is(x, "select.explore")) {

      if (is(x, "summary")) {
        print_summary_select_explore(x, ...)


      } else {
        print_select_explore(x, ...)
      }

    }



    if (is(x, "select.ggm_compare_bf")) {
      print_select_ggm_compare_bf(x, ...)
    }


    if (is(x, "ggm_compare_explore")) {
      if (is(x, "summary.ggm_compare_explore")) {
        print_summary_ggm_compare_bf(x)

      } else {
        print_ggm_compare_bf(x, ...)

      }
    }
  } # end of explore

  if (is(x, "confirm")) {
    if (is(x, "ggm_compare_confirm")) {
      print_ggm_confirm(x , ...)


    } else {
      print_confirm(x, ...)
    }
  }  # end confirm

  # coefficients
  if(is(x, "coef")){

    if(is(x, "summary.coef")){

      print_summary_coef(x, ...)

      } else {

        print_coef(x,...)

      }
  }

if(is(x, "roll_your_own")){

  print_roll_your_own(x, ...)
 }

}
