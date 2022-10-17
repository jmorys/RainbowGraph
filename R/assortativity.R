

mix_mat <- function(atrr_list){
  # atrr_list columns : c1, c2, weight
  mmax <- max(atrr_list[,1:2])
  fillingframe <- dplyr::tibble(c1 = 1:mmax, c2 = 1:mmax, weight = rep(0, mmax))
  atrr_list <- dplyr::bind_rows(atrr_list, fillingframe)
  sum_list <- atrr_list %>% dplyr::group_by(c1, c2) %>% dplyr::summarise(sum = sum(weight), .groups = "drop" )
  mymat <- tidyr::spread(sum_list, c1, sum)
  s <- as.matrix(mymat[,-1])
  s[is.na(s)] <- 0
  s <- s + t(s)
  s <- s/sum(s)
  s
}



assortativity_weighted <- function(graph, val) {
  val <- as.integer(val)
  c1 <- val[igraph::get.edgelist(graph, names = F)[,1]]
  c2 <- val[igraph::get.edgelist(graph, names = F)[,2]]

  atrr_list <- dplyr::tibble(c1 = c1, c2 = c2,
                             weight = igraph::get.edge.attribute(graph, name = "weights"))
  s <- mix_mat(atrr_list)
  sum_s2 <- sum(s %*% s)
  if (dim(s)[1] == 1 | sum_s2 == 1) return(1)
  (sum(diag(s)) - sum_s2) / (1 - sum_s2)
}




#' Calculate local assortativity in parallel
#'
#' @param graph
#' @param val
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
assortativity_local_par <- function(graph, val, alpha = 0.2){
  cl = parallel::makeCluster(parallel::detectCores(logical = FALSE)-1)
  on.exit(parallel::stopCluster(cl))
  graph_l <- igraph::vcount(graph)
  proceed <- graph_l == length(val)
  if (!proceed) stop("length of data differs")
  val <- as.integer(val)
  e2 <- igraph::get.edgelist(graph, names = F)
  e1 <- e2[,1]
  e2 <- e2[,2]
  atrr_list <- dplyr::tibble(c1 = val[e1], c2 = val[e2], weight = rep(1, length(val[e1])) ) # NO WEIGHTS!!!
  s <- mix_mat(atrr_list)
  agg <- sum(s %*% s)
  if (agg == 1) agg <- 0.9999999999
  #sameness bool
  c_bool <- val[c(e1,e2)] == val[c(e2,e1)]
  deg <- as.numeric(igraph::degree(graph)[c(e1,e2)])
  assorts <- rep(0, graph_l)
  cc <- c(e1,e2)
  parallel::clusterExport(cl, c("graph_l", "deg", "graph", "alpha", "cc"), envir = environment())
  assorts = parallel::parSapply(cl, X = 1:graph_l, FUN = function(x){
    pers <- rep(0, graph_l)
    pers[x] <- 1
    Vpage <- igraph::page_rank(graph, personalized = pers, damping =  alpha, weights = NULL)$vector
    Vpage <- Vpage[cc]/deg
    Vpage[Vpage < 0] <- 0
    sum(Vpage[c_bool])/sum(Vpage)
  })

  assorts <- (assorts - agg)/(1 - agg)
  assorts
}
assortativity_local_par <- compiler::cmpfun(assortativity_local_par)



fast_assorts_p_1 <- function(graph, alpha = 0.4){
  cl = parallel::makeCluster(parallel::detectCores(logical = FALSE)-1)
  on.exit(parallel::stopCluster(cl))
  graph_l <- igraph::vcount(graph)
  e2 <- igraph::get.edgelist(graph, names = F)
  e1 <- e2[,1]
  e2 <- e2[,2]
  cc <- c(e1, e2)
  deg <- as.numeric(igraph::degree(graph)[c(e1,e2)])
  parallel::clusterExport(cl, c("graph_l", "deg", "graph", "alpha", "cc"), envir = environment())
  assorts_ori = parallel::parSapply(cl, X = 1:graph_l, FUN = function(x){
    pers <- rep(0, graph_l)
    pers[x] <- 1
    Vpage <- igraph::page_rank(graph, personalized = pers, damping =  alpha, weights = NULL)$vector
    Vpage <- Vpage[cc]/deg
    Vpage
  })
  assorts_ori[assorts_ori < 0] <- 0
  return(list(assorts_ori, e1, e2))
}
# fast_assorts_p_1 <- compiler::cmpfun(fast_assorts_p_1)



fast_assorts_p_2 <- function(val, assorts, e1, e2){
  val <- as.integer(val)
  atrr_list <- dplyr::tibble(c1 = val[e1], c2 = val[e2], weight = rep(1, length(val[e1])) ) # NO WEIGHTS!!!
  s <- mix_mat(atrr_list)
  agg <- sum(s %*% s)
  if (agg == 1) agg <- 0.9999
  #sameness bool
  c_bool <- val[c(e1,e2)] == val[c(e2,e1)]
  c_bool <- c_bool*1
  assorts <- colSums(assorts*c_bool)/colSums(assorts)
  assorts <- (assorts - agg)/(1 - agg)
  assorts
}
# fast_assorts_p_2 <- compiler::cmpfun(fast_assorts_p_2)
