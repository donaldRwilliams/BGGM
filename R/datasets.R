#' post-traumatic stress disorder dataset
#'
#' A dataset containing items that measure Post-traumatic stress disorder symptoms. There are 20 variables (p) and
#' 221 observations (n)
#'
#' \itemize{
#'   \item Intrusive Thoughts
#'   \item Nightmares
#'   \item Flashbacks
#'   \item Emotional cue reactivity
#'   \item Psychological cue reactivity
#'   \item Avoidance of thoughts
#'   \item Avoidance of reminders
#'   \item Trauma-related amnesia
#'   \item Negative beliefs
#'   \item Negative trauma-related emotions
#'   \item Loss of interest
#'   \item Detachment
#'   \item Restricted affect
#'   \item Irritability/anger
#'   \item Self-destructive/reckless behavior
#'   \item Hypervigilance
#'   \item Exaggerated startle response
#'   \item Difficulty concentrating
#'   \item Sleep disturbance
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ptsd
#' @usage data("ptsd")
#' @references
#' Epskamp, S., & Fried, E. I. (2018). A tutorial on regularized partial correlation networks. Psychological methods.
#'  @format A data frame with 221 rows and 20 variables
NULL



#' post-traumatic stress disorder correlation matrix (# 1)
#'
#' A dataset containing items that measure Post-traumatic stress disorder symptoms. There are 16 variables in total.
#' The correlation matrix was estimated from 526 individuals.
#'
#' \itemize{
#'   \item Intrusive Thoughts
#'   \item Nightmares
#'   \item Flashbacks
#'   \item Physiological/psychological reactivity
#'   \item Avoidance of thoughts
#'   \item Avoidance of situations
#'   \item Amnesia
#'   \item Disinterest in activities
#'   \item Feeling detached
#'   \item Emotional numbing
#'   \item Foreshortened future
#'   \item Sleep problems
#'   \item Irritability
#'   \item Concentration problems
#'   \item Hypervigilance
#'   \item Startle response
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ptsd_cor1
#' @examples
#' data(ptsd_cor1)
#' Y <- MASS::mvrnorm(n = 965, mu = rep(0, 16),
#'                    Sigma = ptsd_cor1,  empirical = TRUE)
#' @references
#' Fried, E. I., Eidhof, M. B., Palic, S., Costantini, G., Huisman-van Dijk, H. M., Bockting, C. L., ... & Karstoft, K. I. (2018).
#' Replicability and generalizability of posttraumatic stress disorder (PTSD) networks: a cross-cultural multisite study of PTSD
#' symptoms in four trauma patient samples. Clinical Psychological Science, 6(3), 335-351.
#' @format A (correlation) matrix with 16 variables
NULL

#' post-traumatic stress disorder correlation matrix (# 2)
#'
#' A dataset containing items that measure Post-traumatic stress disorder symptoms. There are 16 variables in total.
#' The correlation matrix was estimated from 365 individuals.
#'
#' \itemize{
#'   \item Intrusive Thoughts
#'   \item Nightmares
#'   \item Flashbacks
#'   \item Physiological/psychological reactivity
#'   \item Avoidance of thoughts
#'   \item Avoidance of situations
#'   \item Amnesia
#'   \item Disinterest in activities
#'   \item Feeling detached
#'   \item Emotional numbing
#'   \item Foreshortened future
#'   \item Sleep problems
#'   \item Irritability
#'   \item Concentration problems
#'   \item Hypervigilance
#'   \item Startle response
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ptsd_cor2
#' @examples
#' data(ptsd_cor2)
#' Y <- MASS::mvrnorm(n = 365, mu = rep(0, 16),
#'                    Sigma = ptsd_cor2,  empirical = TRUE)
#' @references
#' Fried, E. I., Eidhof, M. B., Palic, S., Costantini, G., Huisman-van Dijk, H. M., Bockting, C. L., ... & Karstoft, K. I. (2018).
#' Replicability and generalizability of posttraumatic stress disorder (PTSD) networks: a cross-cultural multisite study of PTSD
#' symptoms in four trauma patient samples. Clinical Psychological Science, 6(3), 335-351.
#' @format A (correlation) matrix with 16 variables
NULL

#' post-traumatic stress disorder correlation matrix (# 3)
#'
#' A dataset containing items that measure Post-traumatic stress disorder symptoms. There are 16 variables in total.
#' The correlation matrix was estimated from 926 individuals.
#'
#' \itemize{
#'   \item Intrusive Thoughts
#'   \item Nightmares
#'   \item Flashbacks
#'   \item Physiological/psychological reactivity
#'   \item Avoidance of thoughts
#'   \item Avoidance of situations
#'   \item Amnesia
#'   \item Disinterest in activities
#'   \item Feeling detached
#'   \item Emotional numbing
#'   \item Foreshortened future
#'   \item Sleep problems
#'   \item Irritability
#'   \item Concentration problems
#'   \item Hypervigilance
#'   \item Startle response
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ptsd_cor3
#' @examples
#' data(ptsd_cor3)
#' Y <- MASS::mvrnorm(n = 926, mu = rep(0, 16),
#'                    Sigma = ptsd_cor3,  empirical = TRUE)
#' @references
#' Fried, E. I., Eidhof, M. B., Palic, S., Costantini, G., Huisman-van Dijk, H. M., Bockting, C. L., ... & Karstoft, K. I. (2018).
#' Replicability and generalizability of posttraumatic stress disorder (PTSD) networks: a cross-cultural multisite study of PTSD
#' symptoms in four trauma patient samples. Clinical Psychological Science, 6(3), 335-351.
#' @format A (correlation) matrix with 16 variables
NULL

#' post-traumatic stress disorder correlation matrix (# 4)
#'
#' A dataset containing items that measure Post-traumatic stress disorder symptoms. There are 16 variables in total.
#' The correlation matrix was estimated from 965 individuals.
#'
#' \itemize{
#'   \item Intrusive Thoughts
#'   \item Nightmares
#'   \item Flashbacks
#'   \item Physiological/psychological reactivity
#'   \item Avoidance of thoughts
#'   \item Avoidance of situations
#'   \item Amnesia
#'   \item Disinterest in activities
#'   \item Feeling detached
#'   \item Emotional numbing
#'   \item Foreshortened future
#'   \item Sleep problems
#'   \item Irritability
#'   \item Concentration problems
#'   \item Hypervigilance
#'   \item Startle response
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ptsd_cor4
#' @examples
#' data(ptsd_cor4)
#' Y <- MASS::mvrnorm(n = 965, mu = rep(0, 16),
#'                    Sigma = ptsd_cor4,  empirical = TRUE)
#' @references
#' Fried, E. I., Eidhof, M. B., Palic, S., Costantini, G., Huisman-van Dijk, H. M., Bockting, C. L., ... & Karstoft, K. I. (2018).
#' Replicability and generalizability of posttraumatic stress disorder (PTSD) networks: a cross-cultural multisite study of PTSD
#' symptoms in four trauma patient samples. Clinical Psychological Science, 6(3), 335-351.
#' @format A (correlation) matrix with 16 variables
NULL

#' 25 Personality items representing 5 factors
#'
#' This data set and the documentation was taken from the \strong{psych} package. Further details can be found in the documentation
#' of the \strong{psych} package.
#'
#' \itemize{
#'   \item \code{A1} Am indifferent to the feelings of others. (q_146)
#'   \item \code{A2} Inquire about others' well-being. (q_1162)
#'   \item \code{A3} Know how to comfort others. (q_1206)
#'   \item \code{A4} Love children. (q_1364)
#'   \item \code{A5} Make people feel at ease. (q_1419)
#'   \item \code{C1} Am exacting in my work. (q_124)
#'   \item \code{C2} Continue until everything is perfect. (q_530)
#'   \item \code{C3} Do things according to a plan. (q_619)
#'   \item \code{C4} Do things in a half-way manner. (q_626)
#'   \item \code{C5} Waste my time. (q_1949)
#'   \item \code{E1} Don't talk a lot. (q_712)
#'   \item \code{E2} Find it difficult to approach others. (q_901)
#'   \item \code{E3} Know how to captivate people. (q_1205)
#'   \item \code{E4} Make friends easily. (q_1410)
#'   \item \code{E5} Take charge. (q_1768)
#'   \item \code{N1} Get angry easily. (q_952)
#'   \item \code{N2} Get irritated easily. (q_974)
#'   \item \code{N3} Have frequent mood swings. (q_1099)
#'   \item \code{N4} Often feel blue. (q_1479)
#'   \item \code{N5} Panic easily. (q_1505)
#'   \item \code{o1} Am full of ideas. (q_128)
#'   \item \code{o2} Avoid difficult reading material.(q_316)
#'   \item \code{o3} Carry the conversation to a higher level. (q_492)
#'   \item \code{o4} Spend time reflecting on things. (q_1738)
#'   \item \code{o5} Will not probe deeply into a subject. (q_1964)
#'   \item \code{gender} Males = 1, Females =2
#'   \item \code{education} 1 = HS, 2 = finished HS, 3 = some college, 4 = college graduate 5 = graduate degree
#'   \item \code{age} age in years
#' }
#'
#' @docType data
#' @keywords datasets
#' @name bfi
#' @references
#' Revelle, W. (2018) psych: Procedures for Personality and Psychological Research, Northwestern University,
#' Evanston, Illinois, USA, https://CRAN.R-project.org/package=psych Version = 1.8.12.
#' @format A data frame with 25 variables and 2800 observations (but with missing data)
NULL
