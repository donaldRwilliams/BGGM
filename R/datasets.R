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



#' Contingencies of Self-Worth Scale (CSWS)
#'
#' A dataset containing items from the Contingencies of Self-Worth Scale (CSWS) scale. There are 35 variables  and
#' 680 observations
#'
#' \itemize{
#'   \item \code{1} When I think I look attractive, I feel good about myself
#'   \item \code{2} My self-worth is based on God’s love
#'   \item \code{3} I feel worthwhile when I perform better than others on a task or skill.
#'   \item \code{4} My self-esteem is unrelated to how I feel about the way my body looks.
#'   \item \code{5} Doing something I know is wrong makes me lose my self-respect
#'   \item \code{6} I don’t care if other people have a negative opinion about me.
#'   \item \code{7} Knowing that my family members love me makes me feel good about myself.
#'   \item \code{8} I feel worthwhile when I have God’s love.
#'   \item \code{9} I can’t respect myself if others don’t respect me.
#'   \item \code{10} My self-worth is not influenced by the quality of my relationships with my family members.
#'   \item \code{11} Whenever I follow my moral principles, my sense of self-respect gets a boost.
#'   \item \code{12} Knowing that I am better than others on a task raises my self-esteem.
#'   \item \code{13} My opinion about myself isn’t tied to how well I do in school.
#'   \item \code{14} I couldn’t respect myself if I didn’t live up to a moral code.
#'   \item \code{15} I don’t care what other people think of me.
#'   \item \code{16} When my family members are proud of me, my sense of self-worth increases.
#'   \item \code{17} My self-esteem is influenced by how attractive I think my face or facial features are.
#'   \item \code{18} My self-esteem would suffer if I didn’t have God’s love.
#'   \item \code{19} Doing well in school gives me a sense of selfrespect.
#'   \item \code{20} Doing better than others gives me a sense of self-respect.
#'   \item \code{21} My sense of self-worth suffers whenever I think I don’t look good.
#'   \item \code{22} I feel better about myself when I know I’m doing well academically.
#'   \item \code{23} What others think of me has no effect on what I think about myself.
#'   \item \code{24} When I don’t feel loved by my family, my selfesteem goes down.
#'   \item \code{25} My self-worth is affected by how well I do when I am competing with others.
#'   \item \code{26} My self-esteem goes up when I feel that God loves me.
#'   \item \code{27} My self-esteem is influenced by my academic performance.
#'   \item \code{28} My self-esteem would suffer if I did something unethical.
#'   \item \code{29} It is important to my self-respect that I have a family that cares about me.
#'   \item \code{30} My self-esteem does not depend on whether or not I feel attractive.
#'   \item \code{31} When I think that I’m disobeying God, I feel bad about myself.
#'   \item \code{32} My self-worth is influenced by how well I do on competitive tasks.
#'   \item \code{33} I feel bad about myself whenever my academic performance is lacking.
#'   \item \code{34} My self-esteem depends on whether or not I follow my moral/ethical principles.
#'   \item \code{35} My self-esteem depends on the opinions others hold of me.
#'   \item \code{gender} "M" (male) or "F" (female)
#'
#' }
#'
#' @note There are seven domains
#'
#' FAMILY SUPPORT: items 7, 10, 16, 24, and 29.
#'
#' COMPETITION: items 3, 12, 20, 25, and 32.
#'
#' APPEARANCE: items 1, 4, 17, 21, and 30.
#'
#' GOD'S LOVE: items 2, 8, 18, 26, and 31.
#'
#' ACADEMIC COMPETENCE: items 13, 19, 22, 27, and 33.
#'
#' VIRTUE: items 5, 11, 14, 28, and 34.
#'
#' APPROVAL FROM OTHERS: items: 6, 9, 15, 23, and 35.
#' @docType data
#' @keywords datasets
#' @name csws
#' @usage data("csws")
#' @examples
#' data("csws")
#'
#' # labels
#' csws_lables <- BGGM:::csws_labels
#'
#' @references
#' Briganti, G., Fried, E. I., & Linkowski, P. (2019). Network analysis of Contingencies of Self-Worth
#' Scale in 680 university students. Psychiatry research, 272, 252-257.
#' @format A data frame with 35 variables and 680 observations (7 point Likert scale)
NULL




#' Toronto Alexithymia Scale (TAS)
#'
#' A dataset containing items from the Toronto Alexithymia Scale (TAS). There are 20 variables  and
#' 1925 observations
#'
#' \itemize{
#'   \item \code{1} I am often confused about what emotion I am feeling
#'   \item \code{2}  It is difficult for me to find the right words for my feelings
#'   \item \code{3} I have physical sensations that even doctors don’t understand
#'   \item \code{4} I am able to describe my feelings easily
#'   \item \code{5} I prefer to analyze problems rather than just describe them
#'   \item \code{6} When I am upset, I don’t know if I am sad, frightened, or angry
#'   \item \code{7}  I am often puzzled by sensations in my body
#'   \item \code{8} I prefer just to let things happen rather than to understand why they turned out that way
#'   \item \code{9}  I have feelings that I can’t quite identify
#'   \item \code{10} Being in touch with emotions is essential
#'   \item \code{11}  I find it hard to describe how I feel about people
#'   \item \code{12} People tell me to describe my feelings more
#'   \item \code{13} I don’t know what’s going on inside me
#'   \item \code{14} I often don’t know why I am angry
#'   \item \code{15} I prefer talking to people about their daily activities rather than their feelings
#'   \item \code{16}  I prefer to watch “light” entertainment shows rather than psychological dramas
#'   \item \code{17} It is difficult for me to reveal my innermost feelings, even to close friends
#'   \item \code{18}  I can feel close to someone, even in moments of silence
#'   \item \code{19}  I find examination of my feelings useful in solving personal problems
#'   \item \code{20} Looking for hidden meanings in movies or plays distracts from their enjoyment
#'   \item \code{gender} "M" (male) or "F" (female)
#'
#' }
#'
#' @note There are three domains
#'
#' Difficulty identifying feelings: items 1, 3, 6, 7, 9, 13, 14
#'
#' Difficulty describing feelings: items 2, 4, 11, 12, 17
#'
#' Externally oriented thinking: items 10, 15, 16, 18, 19
#'
#' @docType data
#' @keywords datasets
#' @name tas
#' @usage data("tas")
#' @examples
#' data("tas")
#'
#' # labels
#' tas_labels <- BGGM:::tas_labels
#' @references
#' Briganti, G., & Linkowski, P. (2019). Network approach to items and domains from
#' the Toronto Alexithymia Scale. Psychological reports.
#' @format A data frame with 20 variables and 1925 observations (5 point Likert scale)
NULL





#' Interpersonal Reactivity Index (IRI)
#'
#' A dataset containing items from the Interpersonal Reactivity Index (IRI; an empathy measure). There are 28 variables  and
#' 1973 observations
#'
#' \itemize{
#'   \item \code{1} I daydream and fantasize, with some regularity, about things that might happen to me.
#'   \item \code{2}  I often have tender, concerned feelings for people less fortunate than me.
#'   \item \code{3} I sometimes find it difficult to see things from the "other guy's" point of view.
#'   \item \code{4} Sometimes I don't feel very sorry for other people when they are having problems.
#'   \item \code{5}  I really get involved with the feelings of the characters in a novel.
#'   \item \code{6} In emergency situations, I feel apprehensive and ill-at-ease.
#'   \item \code{7}  I am usually objective when I watch a movie or play, and I don't often get completely caught up in it.
#'   \item \code{8} I try to look at everybody's side of a disagreement before I make a decision.
#'   \item \code{9}  When I see someone being taken advantage of, I feel kind of protective towards them.
#'   \item \code{10} I sometimes feel helpless when I am in the middle of a very emotional situation.
#'   \item \code{11}  I sometimes try to understand my friends better
#'   by imagining how things look from their perspective
#'   \item \code{12} Becoming extremely involved in a good book or movie is somewhat rare for me.
#'   \item \code{13} When I see someone get hurt, I tend to remain calm.
#'   \item \code{14} Other people's misfortunes do not usually disturb me a great deal.
#'   \item \code{15}  If I'm sure I'm right about something, I don't waste much
#'   time listening to other people's arguments.
#'   \item \code{16}  After seeing a play or movie, I have felt as though I were one of the characters.
#'   \item \code{17}  Being in a tense emotional situation scares me.
#'   \item \code{18}  When I see someone being treated unfairly,
#'   I sometimes don't feel very much pity for them.
#'   \item \code{19} I am usually pretty effective in dealing with emergencies.
#'   \item \code{20} I am often quite touched by things that I see happen.
#'   \item \code{21} I believe that there are two sides to every question and try to look at them both.
#'   \item \code{22} I would describe myself as a pretty soft-hearted person.
#'   \item \code{23} When I watch a good movie, I can very easily put myself in
#'   the place of a leading character
#'   \item \code{24}  I tend to lose control during emergencies.
#'   \item \code{25} When I'm upset at someone, I usually try to "put myself in his shoes" for a while.
#'   \item \code{26} When I am reading an interesting story or novel, I imagine how I would feel if the
#'   events in the story were happening to me.
#'   \item \code{27}  When I see someone who badly needs help in an emergency, I go to pieces.
#'   \item \code{28} Before criticizing somebody, I try to imagine how I would feel if I were in their place.
#'   \item \code{gender} "M" (male) or "F" (female)
#'
#' }
#'
#' @note There are four domains
#'
#' Fantasy: items 1, 5, 7, 12, 16, 23, 26
#'
#' Perspective taking: items 3, 8, 11, 15, 21, 25, 28
#'
#' Empathic concern: items 2, 4, 9, 14, 18, 20, 22
#'
#' Personal distress: items 6, 10, 13, 17, 19, 24, 27,
#'
#' @docType data
#' @keywords datasets
#' @name iri
#' @usage data("iri")
#' @examples
#' data("iri")
#'
#' # labels
#' iri_labels <- BGGM:::iri_labels
#' @references
#'Briganti, G., Kempenaers, C., Braun, S., Fried, E. I., & Linkowski, P. (2018). Network analysis of
#'empathy items from the interpersonal reactivity index in 1973
#'young adults. Psychiatry research, 265, 87-92.
#' @format A data frame with 28 variables and 1973 observations (5 point Likert scale)
NULL




#' Resilience Scale of Adults (RSA)
#'
#' A dataset containing items from the Resilience Scale of Adults (RSA). There are 33 items  and
#' 675 observations
#'
#' \itemize{
#'   \item \code{1}  My plans for the future are
#'   \item \code{2}  When something unforeseen happens
#'   \item \code{3}  My family understanding of what is important in life is
#'   \item \code{4}  I feel that my future looks
#'   \item \code{5}  My goals
#'   \item \code{6}  I can discuss personal issues with
#'   \item \code{7}  I feel
#'   \item \code{8}  I enjoy being
#'   \item \code{9}  Those who are good at encouraging are
#'   \item \code{10} The bonds among my friends
#'   \item \code{11} My personal problems
#'   \item \code{12} When a family member experiences a crisis/emergency
#'   \item \code{13} My family is characterised by
#'   \item \code{14} To be flexible in social settings
#'   \item \code{15} I get support from
#'   \item \code{16} In difficult periods my family
#'   \item \code{17} My judgements and decisions
#'   \item \code{18} New friendships are something
#'   \item \code{19} When needed, I have
#'   \item \code{20} I am at my best when I
#'   \item \code{21} Meeting new people is
#'   \item \code{22} When I am with others
#'   \item \code{23} When I start on new things/projects
#'   \item \code{24} Facing other people, our family acts
#'   \item \code{25} Belief in myself
#'   \item \code{26} For me, thinking of good topics of conversation is
#'   \item \code{27} My close friends/family members
#'   \item \code{28} I am good at
#'   \item \code{29} In my family, we like to
#'   \item \code{30} Rules and regular routines
#'   \item \code{31} In difficult periods I have a tendency to
#'   \item \code{32} My goals for the future are
#'   \item \code{33} Events in my life that I cannot influence
#'   \item \code{gender} "M" (male) or "F" (female)
#'
#' }
#'
#' @note There are 6 domains
#'
#' Planned future: items 1, 4, 5, 32
#'
#' Perception of self: items 2, 11, 17, 25, 31, 33
#'
#' Family cohesion: items 3, 7, 13, 16, 24, 29
#'
#' Social resources: items 6, 9, 10, 12, 15, 19, 27
#'
#' Social Competence: items 8, 14, 18, 21, 22, 26,
#'
#' Structured style: items 23, 28, 30
#'
#' @docType data
#' @keywords datasets
#' @name rsa
#' @usage data("rsa")
#' @examples
#' data("rsa")
#'
#' # labels
#' rsa_labels <- BGGM:::rsa_labels
#' @references
#' Briganti, G., & Linkowski, P. (2019). Item and domain network structures of the Resilience
#' Scale for Adults in 675 university students. Epidemiology and psychiatric sciences, 1-9.
#' @format A data frame with 28 variables and 1973 observations (5 point Likert scale)
NULL

