
#' Simulate regeneration with each cell
#'
#' @param graph undirected graph with edge weight parameter
#' @param survivor_fraction fraction of cells that survive
#' @param expanding_fraction fraction of cells among survivors that multiply and drive regeneration
#'
#' @return list containing vectors of surviving cells, multiplying cells, cells after expansion without non_multiplying cells and finally with non_multiplying cells.
#' The results of this function reference unique ids of cells, not colors.
#'
#' @export
#'
#' @examples
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
simulate_expansion_of_fraction <- function(graph, survivor_fraction=0.3, expanding_fraction=0.1){
  nlen <- igraph::vcount(graph)
  surv_nodes <- sample(1:nlen, ceiling(nlen*survivor_fraction))
  survivors <- rep(NA, nlen)
  survivors[surv_nodes] <- surv_nodes
  ori_nodes <- sample(surv_nodes, ceiling(length(surv_nodes)*expanding_fraction))
  oris <- rep(NA, nlen)
  oris[ori_nodes] <- ori_nodes
  oris <- as.factor(oris)
  cl = parallel::makeCluster(7)
  on.exit(parallel::stopCluster(cl))
  parallel::clusterExport(cl, c("nlen", "graph"), envir = environment())


  non_exp <- surv_nodes[!surv_nodes %in% ori_nodes]
  igraph::E(graph)[incident(graph, non_exp)]$weight <- igraph::E(graph)[incident(graph, non_exp)]$weight * 0.5


  conn = parSapply(cl, X = ori_nodes, FUN = function(x){
    pers <- rep(0, nlen)
    pers[x] <- 1
    igraph::page_rank(graph, personalized = pers, damping =  0.5, weights = NULL)$vector
  })


  expanded <- parallel::parApply(cl, conn,1, function(x){
    order(x, decreasing = T)[1]
  })


  dead <- colMeans(distances(graph, v = sort(ori_nodes)) == Inf) == 1
  expanded[dead] <- NA
  final <- expanded
  final[non_exp] <- max(expanded, na.rm = T) + (1:length(non_exp))
  list(survivors = survivors, origins = oris, expanded = expanded, final = final)
}





#' Simulate regeneration with colors
#'
#' @param graph undirected graph with edge weight parameter
#' @param survivor_fraction fraction of cells that survive
#' @param expanding_fraction fraction of cells among survivors that multiply and drive regeneration
#' @param cols vectors of cell colors. Suggested are of length equal to graph's node count
#' @param with_push whenever to apply repulsing force to multiplying cells
#' @return list containing vectors of surviving cells, multiplying cells, cells after expansion without non_multiplying cells and finally with non_multiplying cells.
#' The results of this function reference colors, with the exception of surviving cells vector. It will also reference cells before
#'
#' @export
#'
#' @examples
simulate_expansion_of_fraction_cheap <- function(graph, survivor_fraction=0.3, expanding_fraction=0.1,
                                                 cols, with_push = F){
  nlen <- igraph::vcount(graph)
  surv_nodes <- sample(1:nlen, ceiling(nlen*survivor_fraction))
  survivors <- rep(NA, nlen)
  survivors[surv_nodes] <- surv_nodes
  ori_nodes <- sample(surv_nodes, ceiling(length(surv_nodes)*expanding_fraction))
  oris <- rep(NA, nlen)
  oris[ori_nodes] <- ori_nodes
  oris <- cols[oris]

  non_exp <- surv_nodes[!surv_nodes %in% ori_nodes]


  conn = sapply(X = sort(na.exclude(unique(oris))), FUN = function(x){
    pers <- rep(0, nlen)
    pers[oris == x] <- 1
    igraph::page_rank(graph, personalized = pers, damping =  0.5, weights = NULL)$vector
  })


  expanded <- apply(conn, 1, function(x){
    order(x, decreasing = T)[1]
  })

  if (with_push) {
    #create directed edge list
    originz <- !is.na(oris)
    out <- igraph::page.rank(graph, personalized = originz*1., directed = F, damping = 0.5)
    dif1_v <- out$vector
    survivors <- !is.na(survivors)
    non_expanding <- (survivors & !originz)
    edge_list <- igraph::get.edgelist(graph, names = F)
    e1s <- dif1_v[edge_list[,1]]
    e2s <- dif1_v[edge_list[,2]]
    edge_list[e1s < e2s, ] <- edge_list[e1s < e2s,2:1 ]

    # prepare edgelist for pebble motion
    edge_params <- tibble::tibble(from = edge_list[,1], to = edge_list[,2], potential = dif1_v[ edge_list[,2]])
    edge_params <- edge_params %>% dplyr::filter(from %in% which(!originz) & to %in% which(!originz))
    edge_params <- edge_params %>% dplyr::arrange(potential)

    non_expanding_update <- non_expanding
    continue <- T
    while (continue) {
      arranged_sub <- edge_params %>% dplyr::filter(from %in% which(non_expanding_update) & to %in% which(!non_expanding_update) )
      substitution <- arranged_sub %>% dplyr::group_by(from) %>% dplyr::filter(row_number() == 1) %>%
        dplyr::group_by(to) %>% dplyr::filter(row_number() == 1)
      non_expanding_update[substitution$from] <- F
      non_expanding_update[substitution$to] <- T
      continue <- dim(substitution)[1] > 0
    }
    non_exp <- which(non_expanding_update)
    survivors <- non_expanding_update | originz
  }


  dead <- colMeans(distances(graph, v = ori_nodes) == Inf) == 1
  expanded[dead] <- NA
  final <- expanded
  final[non_exp] <- cols[max(expanded, na.rm = T) + (1:length(non_exp))]
  list(survivors = survivors, origins = oris, expanded = expanded, final = final,
       n_surv = length(surv_nodes), n_exp = length(ori_nodes))
}




#' Quickly perform simulations and measure assortativities for a given graph
#'
#' @param graph input undirected graph with weight edge attribute
#' @param simulation_parameters survival and expansion parameters as 2 column dataframe
#' @param alpha damping factor for assortativity measurement
#' @param n_try number of repeats for each entry in simulation parameters
#' @param with_push implement repulsing force in simulation
#'
#' @export
#' @importFrom foreach %dopar% %:%
generate_data_fast <- function(graph, simulation_parameters, alpha = 0.7, n_try = 3, with_push = F){
  fast_outs = list()
  for (i in 1:length(alpha)) {
    fast_outs[[i]] <- fast_assorts_p_1(graph, alpha = alpha[i])
  }

  cl <- parallel::makeCluster(7)
  doParallel::registerDoParallel(cl)
  init_size <- igraph::vcount(graph)
  res <- foreach::foreach(sim_params = iterators::iter(simulation_parameters, by = 'row'),
                          .export = c("mix_mat", "simulate_expansion_of_fraction_cheap",
                                      "fast_assorts_p_2"),
                          .packages = c("igraph", "dplyr")) %:% foreach::foreach(try = 1:n_try) %dopar% {
                            survivor_f <- sim_params[1]
                            expanding_f <- sim_params[2]
                            set.seed(as.integer(234 + 456789*survivor_f + 123456*expanding_f + try*123))
                            cols <- sample(1:3, init_size, replace = T)
                            test_cheap <- simulate_expansion_of_fraction_cheap(graph = graph,
                                                                               survivor_fraction = survivor_f,
                                                                               expanding_fraction = expanding_f,
                                                                               cols = cols, with_push = with_push)
                            n_surv <- test_cheap$n_surv
                            n_exp <- test_cheap$n_exp

                            asses <- list()
                            for (i in 1:length(fast_outs)) {
                              ass <- fast_assorts_p_2(test_cheap$final, fast_outs[[i]][[1]], fast_outs[[i]][[2]], fast_outs[[i]][[3]])
                              asses[[i]] <- ass
                            }

                            list(asses, c(n_surv, n_exp))
                          }
  parallel::stopCluster(cl)
  res <- unlist(recursive = F, res)
  return(res)
}
