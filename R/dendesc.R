#' @title
#' Dendrogram Descriptive Measures
#'
#' @description
#' Descriptive measures for analyzing objects of class
#' \code{"\link[stats]{dendrogram}"}.
#'
#' @param dendro
#'   Object of class \code{"\link[stats]{dendrogram}"} as produced by
#'   \code{\link{linkage}()} or by \code{\link[=dendrogram]{as.dendrogram}()}
#'   applied to the hierarchical trees returned by \code{\link[stats]{hclust}()}
#'   and \code{\link[cluster]{agnes}()}.
#' @param prox
#'   Object of class \code{"\link[stats]{dist}"} containing the proximity data
#'   used to build the dendrogram.
#' @param ultr
#'   Object of class \code{"\link[stats]{dist}"} containing the ultrametric
#'   distances in the dendrogram, sorted in the same order as the proximity data
#'   in \code{prox}.
#'
#' @details
#'
#' This package allows the calculation of several descriptive measures for
#' dendrograms, such as normalized tree balance, cophenetic correlation
#' coefficient, normalized mean absolute error, and space distortion ratio.
#'
#' For each node in a dendrogram, its entropy is calculated using the
#' concept of Shannon's entropy, which gives a maximum entropy of 1 to nodes
#' merging subdendrograms with the same number of leaves. The average entropy
#' for all nodes in a dendrogram is called its tree balance. Normalized 
#' tree balance is computed by the \code{\link{ntb}()} function as the ratio
#' between the tree balance of a dendrogram and the minimum tree balance of any
#' dendrogram with the same number of elements. Perfectly balanced dendrograms
#' have a normalized tree balance equal to 1, while binary dendrograms formed
#' chaining one new element at a time have a normalized tree balance equal to 0.
#'
#' To calculate the cophenetic correlation coefficient, the
#' \code{\link[stats]{cor}()} function in the \pkg{stats} package needs that the
#' matrix of ultrametric distances (also known as cophenetic distances) and the
#' matrix of proximity data used to build the corresponding dendrogram, they
#' both have their rows and columns sorted in the same order. When the
#' \code{\link[stats]{cophenetic}()} function is used with objects of class
#' \code{"\link[stats]{hclust}"}, it returns ultrametric matrices sorted in
#' appropriate order. However, when the \code{\link[stats]{cophenetic}()}
#' function is used with objects of class \code{"\link[stats]{dendrogram}"}, it
#' returns ultrametric matrices sorted in the order of dendrogram leaves. The
#' \code{\link{ultrametric}()} function in this package returns ultrametric
#' matrices in appropriate order to calculate the cophenetic correlation
#' coefficient using the \code{\link[stats]{cor}()} function.
#'
#' The space distortion ratio of a dendrogram is computed by the
#' \code{\link{sdr}()} function as the difference between the maximum and
#' minimum ultrametric distances, divided by the difference between the
#' maximum and minimum original distances used to build the dendrogram. Space
#' dilation occurs when the space distortion ratio is greater than 1.
#'
#' @seealso
#' \code{\link{linkage}()} in this package, \code{\link[stats]{hclust}()} in the
#' \pkg{stats} package, and \code{\link[cluster]{agnes}()} in the \pkg{cluster}
#' package for building hierarchical trees.
#'
#' @examples
#' ## distances between 21 cities in Europe
#' data(eurodist)
#'
#' ## comparison of dendrograms in terms of the following descriptive mesures:
#' ## - normalized tree balance
#' ## - cophenetic correlation coefficient
#' ## - normalized mean absolute error
#' ## - space distortion ratio
#'
#' ## single linkage (call to the mdendro package)
#' dendro1 <- linkage(eurodist, method="single")
#' ntb(dendro1)          # 0.2500664
#' ultr1 <- ultrametric(dendro1)
#' cor(eurodist, ultr1)  # 0.7842797
#' mae(eurodist, ultr1)  # 0.6352011
#' sdr(eurodist, ultr1)  # 0.150663
#'
#' ## complete linkage (call to the stats package)
#' dendro2 <- as.dendrogram(hclust(eurodist, method="complete"))
#' ntb(dendro2)          # 0.8112646
#' ultr2 <- ultrametric(dendro2)
#' cor(eurodist, ultr2)  # 0.735041
#' mae(eurodist, ultr2)  # 0.8469728
#' sdr(eurodist, ultr2)  # 1
#'
#' ## unweighted arithmetic linkage (UPGMA)
#' dendro3 <- linkage(eurodist, method="arithmetic", weighted=FALSE)
#' ntb(dendro3)          # 0.802202
#' ultr3 <- ultrametric(dendro3)
#' cor(eurodist, ultr3)  # 0.7279432
#' mae(eurodist, ultr3)  # 0.294578
#' sdr(eurodist, ultr3)  # 0.5066903
#'
#' ## unweighted geometric linkage
#' dendro4 <- linkage(eurodist, method="geometric", weighted=FALSE)
#' ntb(dendro4)          # 0.7531278
#' ultr4 <- ultrametric(dendro4)
#' cor(eurodist, ultr4)  # 0.7419569
#' mae(eurodist, ultr4)  # 0.2891692
#' sdr(eurodist, ultr4)  # 0.4548112
#'
#' @name dendesc
NULL

#  Normalized Tree Balance  ----------------------------------------------------

#' @describeIn dendesc
#' Returns a number between 0 and 1 representing the normalized tree balance of
#' the input dendrogram.
#'
#' @export
ntb <- function(dendro) {
  tuple <- recbalance(dendro)
  if (tuple$junctions == 0) {
    return(0.0)
  } else {
    b <- tuple$entropy / tuple$junctions
    bmin <- minbalance(members=stats::nobs(dendro))
    return((b - bmin) / (1.0 - bmin))
  }
}

recbalance <- function(dendro) {
  tuple <- list(entropy=0.0, junctions=0)
  if (!stats::is.leaf(dendro)) {
    tuple$entropy <- entropy(dendro)
    tuple$junctions <- 1
    for (i in 1:length(dendro)) {
      subtuple <- recbalance(dendro[[i]])
      tuple$entropy <- tuple$entropy + subtuple$entropy
      tuple$junctions <- tuple$junctions + subtuple$junctions
    }
  }
  return(tuple)
}

entropy <- function(dendro) {
  h <- 0.0
  if (!stats::is.leaf(dendro)) {
    for (i in 1:length(dendro)) {
      p <- stats::nobs(dendro[[i]]) / stats::nobs(dendro)
      h <- h - p * log(p, base=length(dendro))
    }
  }
  return(h)
}

minbalance <- function(members) {
  terms <- vapply(2:(members - 1), FUN=minbalanceterm, FUN.VALUE=double(1))
  return((sum(terms) + log2(members)) / (members - 1))
}

minbalanceterm <- function(n) {
	log2(n) / (n + 1)
}

#  Ultrametric Distances  ------------------------------------------------------

#' @describeIn dendesc
#' Returns an object of class \code{"\link[stats]{dist}"} containing the
#' ultrametric distance matrix sorted in the same order as the proximity matrix
#' used to build the corresponding dendrogram.
#'
#' @export
ultrametric <- function(dendro) {
  coph <- stats::cophenetic(dendro)
  ord <- order(stats::order.dendrogram(dendro))
  ultr <- stats::as.dist(as.matrix(coph)[ord, ord])
  return(ultr)
}

#  Normalized Mean Absolute Error  ---------------------------------------------

#' @describeIn dendesc
#' Returns the normalized mean absolute error.
#'
#' @export
mae <- function(prox, ultr) {
  return(sum(abs(prox - ultr)) / sum(abs(prox)))
}

#  Space Distortion Ratio  -----------------------------------------------------

#' @describeIn dendesc
#' Returns the space distortion ratio.
#'
#' @export
sdr <- function(prox, ultr) {
  return((max(ultr) - min(ultr)) / (max(prox) - min(prox)))
}
