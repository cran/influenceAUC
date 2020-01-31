#' @title  Influence Functions On AUC
#'
#' @description Provide two sample versions (DEIF and SIF) of influence function on the AUC.
#'
#' @param score A vector containing the predictions (continuous scores) assigned by classifiers; Must be numeric.
#' @param binary A vector containing the true class labels 1: positive and 0: negative. Must have the same dimensions as 'score.'
#' @param threshold A numeric value determining the threshold to distinguish influential observations from normal ones; Must lie between 0 and 1; Defaults to 0.5.
#' @param hypothesis Logical which controls the evaluation of SIF under asymptotic distribution.
#' @param testdiff A numeric value determining the difference in the hypothesis testing; Must lie between 0 and 1; Defaults to 0.5.
#' @param alpha A numeric value determining the significance level in the hypothesis testing; Must lie between 0 and 1; Defaults to 0.05.
#' @param name A vector comprising the appellations for observations; Must have the same dimensions as 'score'.
#' @return A list of objects including (1) `output`: a list of results with `AUC` (numeric), `SIF` (a list of dataframes) and `DEIF` (a list of dataframes)); (2) `rdata`: a dataframe of essential results for visualization
#' (3) `threshold`: a used numeric value to distinguish influential observations from normal ones; (4) `test_output`: a list of dataframes for hypothesis testing result; (5) `test_data`: a dataframe of essential results in hypothesis testing for visualization
#' (6) `testdiff`: a used numeric value to determine the difference in the hypothesis testing; (7) `alpha`: a used nuermic value to determine the significance level.
#' @import dplyr ggplot2 ggrepel ROCR
#' @importFrom stats ecdf
#' @importFrom stats qnorm
#' @importFrom methods show
#' @importFrom graphics par
#'
#' @details Apply two sample versions of influence functions on AUC:
#' \itemize{
#'  \item deleted empirical influence function (DEIF)
#'  \item sample influence function (SIF)
#' }
#' The concept of influence function focuses on the deletion diagnostics; nevertheless, such techniques may face masking effect due to multiple influential observations.
#' To thoroughly investigate the potential cases in binary classification, we suggest end-users to apply \code{\link{ICLC}} and \code{\link{LAUC}} as well. For a complete discussion of these functions, please see the reference.
#' @seealso \code{\link{ICLC}}, \code{\link{LAUC}}
#' @export
#' @author  Bo-Shiang Ke and Yuan-chin Ivan Chang
#' @references Ke, B. S., Chiang, A. J., & Chang, Y. C. I. (2018). Influence Analysis for the Area Under the Receiver Operating Characteristic Curve. Journal of biopharmaceutical statistics, 28(4), 722-734.
#' @examples
#' library(ROCR)
#' data("ROCR.simple")
#' # print out IAUC results directly
#' IAUC(ROCR.simple$predictions,ROCR.simple$labels,hypothesis = "True")
#'
#' data(mtcars)
#' glmfit <- glm(vs ~ wt + disp, family = binomial, data = mtcars)
#' prob <- as.vector( predict(glmfit, newdata = mtcars,type = "response"))
#' output <- IAUC(prob, mtcars$vs, threshold = 0.3, testdiff = 0.3,
#'                hypothesis = TRUE, name = rownames(mtcars))
#' # Show results
#' print(output)
#' # Visualize results
#' plot(output)

IAUC <- function(score, binary, threshold = 0.5, hypothesis = FALSE, testdiff = 0.5, alpha = 0.05, name = NULL) {
  # score: the continuous score assigned by the classifier
  # binary: 1 = positive outcome and 0 = negative outcome
  # Include an error if score is not numeric
  if(!is.numeric(score))
    stop("Classification scores must be numeric")

  # Include an error if threshold is not between 0 and 1
  if(!is.numeric(threshold) | threshold > 1 | threshold < 0)
    stop("threshold must be numeric and between 0 and 1")

  # Include an error if testdiff is not between 0 and 1
  if(!is.numeric(testdiff) | testdiff > 1 | testdiff < 0)
    stop("testdiff must be numeric and between 0 and 1")

  # Include an error if alpha is not between 0 and 1
  if(!is.numeric(alpha) | alpha > 1 | alpha < 0)
    stop("alpha must be numeric and between 0 and 1")

  # Include an error if binary is not 0 or 1
  if(!is.numeric(binary) |  sum(!binary %in% c(0,1)) > 0)
    stop("binary response must be 0 or 1")

  # positions of two outcomes
  neg <- which(binary == 0)
  pos <- which(binary == 1)

  # scores the two outcomes
  sneg <- score[neg]
  spos <- score[pos]

  # numbers of two outcomes and overall data
  m <- length(neg)
  n <- length(pos)
  len <- m + n

  # overall AUC via placement value formulas
  # AUC =  mean(F(spos)) =  mean( 1-G(sneg) )
  # F. and S are the cdf and sruvival function of positive (case) scores
  # G is the cdf of negative (cntrol) scores
  G  <- ecdf(sneg)
  F. <- ecdf(spos)
  AUC <- mean(G(spos))

  # pauc: AUCs of each positive value and the whole negative outcome
  pauc <- G(spos) - AUC
  # nauc: AUCs of each negative value and the whole positive outcome
  nauc <- 1 - F.(sneg) - AUC

  #specify the influences vectors for SIF and DEIF
  SIF <- DEIF <- rep(NA, len)

  #SIF
  SIF[c(neg)] <- nauc
  SIF[c(pos)] <- pauc

  #DEIF
  DEIF[c(neg)] <- (m / (m - 1)) * nauc
  DEIF[c(pos)] <- (n / (n - 1)) * pauc

  SIFdetect <- which(abs(SIF) > threshold)
  DEIFdetect <- which(abs(DEIF) > threshold)

  SIFpos <-  intersect(pos, SIFdetect)
  SIFneg <-  intersect(neg, SIFdetect)
  DEIFpos <-  intersect(pos, DEIFdetect)
  DEIFneg <-  intersect(neg, DEIFdetect)

  sif.p <- cbind(Index = SIFpos, influence = SIF[SIFpos])
  sif.n <- cbind(Index = SIFneg, influence = SIF[SIFneg])
  deif.p <- cbind(Index = DEIFpos, influence = DEIF[DEIFpos])
  deif.n <- cbind(Index = DEIFneg, influence = DEIF[DEIFneg])

  # Collect essential elements as a dataframe
  newdata <- data.frame(score, binary, SIF, DEIF) %>%
    mutate(Outcome = ifelse(binary == 1, "Positive", "Negative")) %>%
    mutate(Index = 1:len, sifpoten = "", deifpoten = "")

  if(length(name) != 0) {
    row.names(sif.p) <- name[SIFpos]
    row.names(sif.n) <- name[SIFneg]
    row.names(deif.p) <- name[DEIFpos]
    row.names(deif.n) <- name[DEIFneg]
    row.names(newdata) <- name
  }

  sif_output <- list(Pos = sif.p, Neg = sif.n)
  deif_output <- list(Pos = deif.p, Neg = deif.n)

  newdata$sifpoten[SIFdetect] <- row.names(newdata)[SIFdetect]
  newdata$deifpoten[DEIFdetect] <- row.names(newdata)[DEIFdetect]

  result <- list(output = list(AUC = AUC,
                               SIF = sif_output,
                               DEIF = deif_output),
                 rdata = newdata,
                 threshold = threshold,
                 test_output = NULL,
                 test_data = NULL,
                 testdiff = testdiff,
                 alpha = alpha)

  if(hypothesis) {
    # Hypothesis testing: asyptotic distribution for SIF to define cutting point
    newSIF <- -SIF
    asympvar <- sum(nauc^2)/m^2 + sum(pauc^2)/n^2

    pivot <- (newSIF - testdiff) / sqrt(asympvar)

    testSIFdetect <- which(pivot > qnorm(1 - alpha))
    testSIFpos <-  intersect(pos, testSIFdetect)
    testSIFneg <-  intersect(neg, testSIFdetect)

    testsif.p <- cbind(Index = testSIFpos, pivot.quantity = pivot[testSIFpos])
    testsif.n <- cbind(Index = testSIFneg, pivot.quantity = pivot[testSIFneg])

    if(length(name) != 0) {
      row.names(testsif.p) <- name[testSIFpos]
      row.names(testsif.n) <- name[testSIFneg]
      row.names(newdata) <- name
    }
    testsif_output <- list(Pos = testsif.p, Neg = testsif.n)

    testdata <- data.frame(score, binary, pivot) %>%
      mutate(Outcome = ifelse(binary == 1,"Positive", "Negative")) %>%
      mutate(Index = 1:len, testpoten = "")
    testdata$testpoten[testSIFdetect] <- row.names(testdata)[testSIFdetect]

    result$test_output <- list(Testing = testsif_output)
    result$test_data <- testdata
  }

  class(result) <- "IAUC"

  return(result)
}

#' @title Local Influence Approaches On AUC
#'
#' @description Apply local influence approaches in terms of slope and curvature on the AUC to quantify the impacts of all observations simultaneously.
#'
#' @param score A vector containing the predictions (continuous scores) assigned by classifiers; Must be numeric.
#' @param binary A vector containing the true class labels 1: positive and 0: negative. Must have the same dimensions as 'score.'
#' @param threshold A numeric value determining the threshold to distinguish influential observations from normal ones; Must lie between 0 and 1; Defaults to 0.2.
#' @param name A vector comprising the appellations for observations; Must have the same dimensions as 'score.'
#' @return A list of objects including (1) `output`: a list of results with `AUC` (numeric), `Slope` (a list of dataframes) and `Curvature` (a list of dataframes)); (2) `rdata`: a dataframe of essential results for visualization
#' (3) `threshold`: a used numeric value to distinguish influential observations from normal ones.
#' @details The influence functions on the AUC focus on the deletion diagnostics; however, such approaches may encounter the masking effect. Rather than dealing with single observations
#' once at a time, local influence methods address this issue by finding the weighted direction of all observations accompanied by the greatest (magnitude) slope and curvature. From the explicit formula based on
#' the slope, local influence methods may face the imbalanced data effect. To thoroughly investigate the potential observation in binary classification, we suggest end-users to apply \code{\link{ICLC}} and \code{\link{IAUC}} as well.
#' For a complete discussion of these functions, please see the reference.
#' @import dplyr geigen ggplot2 ggrepel ROCR
#' @importFrom stats ecdf
#' @importFrom methods show
#'
#' @seealso \code{\link{ICLC}}, \code{\link{IAUC}}
#' @export
#' @author  Bo-Shiang Ke and Yuan-chin Ivan Chang
#' @references Ke, B. S., Chiang, A. J., & Chang, Y. C. I. (2018). Influence Analysis for the Area Under the Receiver Operating Characteristic Curve. Journal of biopharmaceutical statistics, 28(4), 722-734.
#' @examples
#' library(ROCR)
#' data("ROCR.simple")
#' # print out LAUC results directly
#' LAUC(ROCR.simple$predictions,ROCR.simple$labels)
#'
#' data(mtcars)
#' glmfit <- glm(vs ~ wt + disp, family = binomial, data = mtcars)
#' prob <- as.vector(predict(glmfit, newdata = mtcars, type = "response"))
#' output <- LAUC(prob, mtcars$vs, name = rownames(mtcars))
#' # Show results
#' print(output)
#' # Visualize results
#' plot(output)
LAUC <- function(score, binary, threshold = 0.2, name = NULL) {
  # score: the continuous score assigned by the classifier
  # binary: 1 = positive outcome and 0 = negative outcome
  # Include an error if score is not numeric
  if(!is.numeric(score))
    stop("Classification scores must be numeric")

  # Include an error if threshold is not between 0 and 1
  if(!is.numeric(threshold) | threshold > 1 | threshold < 0)
    stop("threshold must be numeric and between 0 and 1")

  # Include an error if binary is not 0 or 1
  if(!is.numeric(binary) | sum(!binary %in% c(0, 1)) > 0)
    stop("binary must be 0 or 1")

  # positions of two outcomes
  neg <- which(binary == 0)
  pos <- which(binary == 1)

  # scores the two outcomes
  sneg <- score[neg]
  spos <- score[pos]

  # numbers of two outcomes and overall data
  m <- length(neg)
  n <- length(pos)
  len <- m + n

  # overall AUC via placement value formulas
  # AUC =  mean(F(spos)) =  mean( 1-G(sneg) )
  # F. and S are the cdf and sruvival function of positive (case) scores
  # G is the cdf of negative (cntrol) scores
  G  <- ecdf(sneg)
  F. <- ecdf(spos)
  AUC <- mean(G(spos))

  # pauc: AUCs of each positive value and the whole negative outcome
  pauc <- G(spos) - AUC
  # nauc: AUCs of each negative value and the whole positive outcome
  nauc <- 1 - F.(sneg) - AUC

  # Based on the relation with SIF to obtain
  # the direction with the maximum slope.
  etadot <- rep(NA, len) # slope, equaiton (14)
  etadot[neg] <- c(nauc / m)
  etadot[pos] <- c(pauc / n)
  unitslope <- etadot / sqrt(sum(etadot^2))

  # derive the tau: each negative versus the whole positive
  negcount <- rep(NA, m)
  poscount <- rep(NA, n)

  for(j in 1:m) negcount[j] <- sum(sneg[j] <= spos)
  for(i in 1:n) poscount[i] <- sum(sneg <= spos[i])

  tau <- c(negcount,poscount)
  delta <- c(rep(n, m), rep(m, n))
  capA <- capB <- matrix(0, nrow = len, ncol = len)

  for(j in 1:m) {
    for(i in 1:n) {
      capA[j, m + i] <- 2 * (sneg[j] <= spos [i])
      capB[j, m + i] <- 2
    }
  }

  # E = eta_double_dot
  capE <- (capA - AUC * capB) / (m * n) - (2 * (tau - AUC * delta) %*% t(delta)) / (m^2 * n^2)

  # k = sqrt( 1 + t(etadot) %*% etadot  )
  k <- as.numeric(sqrt(1 + t(etadot) %*% etadot))
  capD <- k * (diag(len) + etadot %*% t(etadot))

  # solve the gerneralized eigenfunction
  ge <- geigen(capE, capD, symmetric = TRUE)

  eigval <- ge$values # relate to the curvature
  eigvec <- ge$vectors # relate to the direction

  # identify the curvature with greatest magnitude
  maxeigval <- which.max(abs(eigval))
  maxvec <- eigvec[, maxeigval]

  curv <- rep(NA, len)
  curv[neg] <- maxvec[1:m]
  curv[pos] <- maxvec[(m + 1):(len)]
  unitcurvature <- curv / sqrt(sum(curv^2))

  slopedetect <- which(abs(unitslope) > threshold)
  curvaturedetect <- which(abs(unitcurvature) > threshold)

  slopepos <- intersect(pos, slopedetect)
  slopeneg <- intersect(neg, slopedetect)
  curvaturepos <- intersect(pos, curvaturedetect)
  curvatureneg <- intersect(neg, curvaturedetect)

  slope.p <- cbind(Index = slopepos, influence = unitslope[slopepos])
  slope.n <- cbind(Index = slopeneg, influence = unitslope[slopeneg])
  curvature.p <- cbind(Index = curvaturepos,influence = unitcurvature[curvaturepos])
  curvature.n <- cbind(Index = curvatureneg,influence = unitcurvature[curvatureneg])

  # Collect essential elements as a dataframe
  newdata <- data.frame(score, binary, unitslope, unitcurvature) %>%
    mutate(Outcome = ifelse(binary == 1, "Positive", "Negative")) %>%
    mutate(Index = 1:len, slopepoten = "", curvaturepoten = "")

  if(length(name) != 0) {
    row.names(slope.p) <- name[slopepos]
    row.names(slope.n) <- name[slopeneg]
    row.names(curvature.p) <- name[curvaturepos]
    row.names(curvature.n) <- name[curvatureneg]
    row.names(newdata) <- name
  }

  slope_output <- list(Pos = slope.p, Neg = slope.n)
  curvature_output <- list(Pos = curvature.p, Neg = curvature.n)

  newdata$slopepoten[slopedetect] <- row.names(newdata)[slopedetect]
  newdata$curvaturepoten[curvaturedetect] <- row.names(newdata)[curvaturedetect]

  result <- list(output = list(AUC = AUC,
                               Slope = slope_output,
                               Curvature = curvature_output),
                 rdata = newdata,
                 threshold = threshold)
  class(result) <- "LAUC"

  return(result)
}

#' @title Cumulative Lift Charts
#'
#' @description  Show the existence and approximate locations of influential observations
#'     in binary classification through modified cumulative lift charts.
#'
#' @param score A vector containing the predictions (continuous scores) assigned by classifiers; Must be numeric.
#' @param binary A vector containing the true class labels 1: positive and 0: negative. Must have the same dimensions as 'score.'
#' @param prop A numeric value determining the proportion; Must lie between 0 and 1; Defaults to 0.2.
#' @return A list of ggplot2 objects
#' @details There are two types of influential cases in binary classification:
#' \itemize{
#'  \item positive cases with relatively lower scores - negative cumulative lift chart (NCLC)
#'  \item negative cases with relatively higher scores - positive cumulative lift chart  (PCLC)
#' }
#' Each cumulative lift chart (PCLC or NCLC) identifies one type of influential observations and mark with red dotted lines. Based on the characteristics
#' of two types of influential cases, identifying them require to search the highest and lowest proportions of 'score.'
#'
#' Graphical approaches only reveal the existence and approximate locations of influential observations; it would be better to include some quantities to measure their impacts
#' to the interested parameter. To fully investigate the potential observation in binary classification, we suggest end-users to apply two quantification
#' methods \code{\link{IAUC}} and  \code{\link{LAUC}} as well. For a complete discussion of these functions, please see the reference.
#' @import ggplot2 dplyr ROCR
#' @importFrom methods show
#'
#' @seealso \code{\link{IAUC}}, \code{\link{LAUC}}
#' @export
#' @author  Bo-Shiang Ke and Yuan-chin Ivan Chang
#' @references Ke, B. S., Chiang, A. J., & Chang, Y. C. I. (2018). Influence Analysis for the Area Under the Receiver Operating Characteristic Curve. Journal of biopharmaceutical statistics, 28(4), 722-734.
#' @examples
#' library(ROCR)
#' data("ROCR.simple")
#' output <- ICLC(ROCR.simple$predictions,ROCR.simple$labels)
#' plot(output)
#' # Customize a text size for NCLC
#' library(ggplot2)
#' output$NCLC + theme(text = element_text(size = 20))
#'
#' data(mtcars)
#' glmfit <- glm(vs ~ wt + disp, family = binomial, data = mtcars)
#' prob <- as.vector(predict(glmfit, newdata = mtcars, type = "response"))
#' plot(ICLC(prob, mtcars$vs, 0.5))
ICLC <- function(score, binary, prop = 0.2) {
  # score: the continuous score assigned by the classifier
  # binary: 1 = positive outcome and 0 = negative outcome
  # prop: display the influential cases within the chosen proportions of the
  #         highest and lowest scores;between 0 and 1. The default is set to 0.2

  # Include an error if score is not numeric
  if(!is.numeric(score))
    stop("Classification scores must be numeric")

  # Include an error if binary is not 0 or 1
  if(!is.numeric(binary) |  sum(!binary %in% c(0, 1)) > 0)
    stop("binary must be 0 or 1")

  # Include an error if prop is not between 0 and 1
  if(!is.numeric(prop) | prop > 1 | prop < 0)
    stop("prop must be numeric and between 0 and 1")

  # Create a positive cumulative lift chart
  pclc <- create_lift_chart(
    score = score,
    binary = binary,
    prop = prop,
    xlab = expression(paste("Positive Predictive Proportion ", beta)),
    ylab = expression(L[r]),
    title = "Positive Cumulative Lift Chart")

  # Create a negative cumulative lift chart
  nclc <- create_lift_chart(
    score = -score,
    binary = ifelse(binary == 1, 0, 1),
    prop = prop,
    xlab = expression(paste("Negative Predictive Proportion ", 1 - beta)),
    ylab = expression(L[l]),
    title = "Negative Cumulative Lift Chart")

  result <- list(PCLC = pclc, NCLC = nclc)
  class(result) <- "ICLC"

  return(result)
}
