# following https://blogs.rstudio.com/ai/posts/2018-11-12-uncertainty_estimates_dropout/




layer_concrete_dropout <- keras::Layer(
  "Layer_Concrete_Dropout",
#' Layer wrapper with concrete dropout
#'
#' @param layer layer around which wrapper is built
#' @param weight_regularizer importance of weight regularizer in final loss - float
#' @param dropout_regularizer importance of dropout regularizer in final loss - float
#' @param init_min minimum value of dropout probability
#' @param init_max maximum value of dropout probability
#' @param is_mc_dropout boolean
#' @param supports_masking boolean
#'
#' @return custom layer
#' @export
  initialize = function(layer,
                        weight_regularizer = 1e-6,
                        dropout_regularizer = 1e-6,
                        init_min = 0.1,
                        init_max = 0.02,
                        is_mc_dropout = T,
                        supports_masking = TRUE
  ) {
    super$initialize()
    self$lay <- layer


    self$weight_regularizer <- weight_regularizer
    self$dropout_regularizer <- dropout_regularizer
    self$is_mc_dropout <- is_mc_dropout
    self$init_min <- keras::k_log(init_min) - keras::k_log(1 - init_min)
    self$init_max <- keras::k_log(init_max) - keras::k_log(1 - init_max)


  },



  build = function(input_shape) {
    self$p_logit <- super$add_weight(
      name = "p_logit",
      shape = keras::shape(1),
      initializer = keras::initializer_random_uniform(self$init_min, self$init_max),
      trainable = TRUE
    )


    if (!self$lay$built) {
      self$lay$build(input_shape)
      self$lay$built <- T
    }

  },

  upp_weights = function(input_shape){
    self$p <- keras::k_sigmoid(self$p_logit)

    input_dim <- input_shape[[2]]


    weight <- self$lay$variables[[1]]
    kernel_regularizer <- self$weight_regularizer *
      keras::k_sum(keras::k_square(weight)) /
      (1 - self$p)


    dropout_regularizer <- self$p * keras::k_log(self$p)
    dropout_regularizer <- dropout_regularizer +
      (1 - self$p) * keras::k_log(1 - self$p)
    dropout_regularizer <- dropout_regularizer *
      self$dropout_regularizer *
      keras::k_cast(input_dim, keras::k_floatx())

    keras::k_sum(kernel_regularizer + dropout_regularizer)
  },


  concrete_dropout = function(x) {
    eps <- keras::k_cast_to_floatx(keras::k_epsilon())
    temp <- 0.1

    unif_noise <- keras::k_random_uniform(shape = keras::k_shape(x))

    drop_prob <- keras::k_log(self$p + eps) -
      keras::k_log(1 - self$p + eps) +
      keras::k_log(unif_noise + eps) -
      keras::k_log(1 - unif_noise + eps)
    drop_prob <- keras::k_sigmoid(drop_prob / temp)
    random_tensor <- 1 - drop_prob

    retain_prob <- 1 - self$p

    x <- x * random_tensor
    x <- x / retain_prob
    x
  },



  call = function(inputs, mask = NULL, training = NULL) {
    reg <- self$upp_weights(inputs$shape)
    self$add_loss(reg)
    self$lay$call(self$concrete_dropout(inputs))
  }
)














#' Calculate Bhattacharyya distance to reference distribution
#'
#' @param ref_dist reference distribution as datapoints
#' @param measure_dists list of distributions
#' @param which_ass which assortativity to read from data
#'
#' @return distance value for each measure distribution
#' @export
#'
get_Bhattacharyya_dist <- function(ref_dist, measure_dists, which_ass = 2){
  d_true <- stats::density(ref_dist[[which_ass]], bw =  0.06, from = -1, to = 1, n = 1024)
  dif_x <- diff(d_true$x)[1]
  d_true <- d_true$y

  bhattach_dist <- purrr::map_dbl(measure_dists, function(x) {
    d_try <- stats::density(x[[1]][[which_ass]],  bw =  0.06, from = -1, to = 1, n = 1024)$y
    BC <- sum(sqrt(d_true * d_try) * dif_x)
    - log(BC)
  })
  return(bhattach_dist)
}








#' Compile assortativities list into df
#'
#' Outputs parameters of local assortativity distributions in the tibble format.
#' Includes moments of diistributions up to kurtosis, quantiles and correlation between distributions.
#'
#' @param assortativities_list list outputed by generate_data function
#'
#' @return Tibble containig distribution parameters of each simulation
#' @export
#'
compile_assortativity_data <- function(assortativities_list){
  r <- (seq(0.1, 0.47^(1/1.8), length.out = 10))^1.8
  q_probes <- sort(c(r,0.5,1 - r))

  #extract patameters of assortativity distributions to estimate their shape
  hists <- purrr::map_dfr(assortativities_list, function(z){
    z <- z[[1]]
    vals <- purrr::map_dfr(z, function(x){
      mean = mean(x)
      var = stats::var(x)
      skew = moments::skewness(x)
      kurt = moments::kurtosis(x)
      quants <- stats::quantile(x, probs = q_probes)
      c1 <- c(mean, var, skew, kurt)
      names(c1) <- c("mean", "var", "skew", "kurt")
      return(c(c1, quants))
    }, .id = "num")
    return(vals)
  }, .id = "assort")

  hists %<>% tidyr::pivot_wider(id_cols = .data$assort, values_from = c(-.data$assort, -.data$num), names_from = .data$num) %>%
    dplyr::select(-.data$assort)

  # corretlation between assortativities obtained with different damping factors
  wh_cors <- purrr::map_dfr(assortativities_list, function(z){
    z <- z[[1]]
    mat <- matrix(unlist(z), nrow = 2)
    mat <- t(mat)
    v_cors <- stats::cor(mat, method = "s")
    vals <- v_cors[upper.tri(v_cors)]
    names(vals) <- 1:length(vals)
    vals
  })
  comb_frame <- dplyr::bind_cols(hists, wh_cors)
  return(comb_frame)
}


as_metrics_df = function(history) {
  # create metrics data frame
  df <- as.data.frame(history$metrics)

  # pad to epochs if necessary
  pad <- history$params$epochs - nrow(df)
  pad_data <- list()
  for (metric in history$params$metrics)
    pad_data[[metric]] <- rep_len(NA, pad)
  df <- rbind(df, pad_data)

  # return df
  df
}











#' Create model
#'
#' @param n_train number of training samples
#' @param input_dim dimensionality of input
#' @param output_dim dimensionality of output
#' @param l prior length-scale
#'
#' @return
#' @export
#'

model_b <-  function(n_train = 3000, input_dim = 20, output_dim = 1, l = 1e-2){

  leaky_relu <- function(x){
    keras::activation_relu(x, alpha = 0.1)
  }

  # initial value for weight regularizer
  wd <- l^2/n_train
  # initial value for dropout regularizer
  dd <- 2/n_train
  # bayesian dropout model
  input <- keras::layer_input(shape = input_dim)
  output <- input %>% layer_concrete_dropout(
    layer = keras::layer_dense(units = 512, activation = leaky_relu ,
                               kernel_constraint = keras::constraint_maxnorm(max_value = 1)),
    weight_regularizer = wd,
    dropout_regularizer = dd
  ) %>% layer_concrete_dropout(
    layer = keras::layer_dense(units = 256, activation = leaky_relu ,
                               kernel_constraint = keras::constraint_maxnorm()),
    weight_regularizer = wd,
    dropout_regularizer = dd
  ) %>% layer_concrete_dropout(
    layer = keras::layer_dense(units = 256, activation = leaky_relu ,
                               kernel_constraint = keras::constraint_maxnorm()),
    weight_regularizer = wd,
    dropout_regularizer = dd
  ) %>% layer_concrete_dropout(
    layer = keras::layer_dense(units = 256, activation = "elu" ,
                               kernel_constraint = keras::constraint_maxnorm()),
    weight_regularizer = wd,
    dropout_regularizer = dd
  )
  mean <- output %>% layer_concrete_dropout(
    layer = keras::layer_dense(units = output_dim ,
                               kernel_constraint = keras::constraint_maxnorm()),
    weight_regularizer = wd,
    dropout_regularizer = dd
  )
  log_var <- output %>% layer_concrete_dropout(
    keras::layer_dense(units = output_dim ,
                       kernel_constraint = keras::constraint_maxnorm()),
    weight_regularizer = wd,
    dropout_regularizer = dd
  )

  output <- keras::layer_concatenate(list(mean, log_var))
  model <- keras::keras_model(input, output)

  heteroscedastic_loss <- function(y_true, y_pred) {
    mean <- y_pred[, 1:output_dim]
    log_var <- y_pred[, (output_dim + 1):(output_dim * 2)]
    precision <- keras::k_exp(-log_var)
    keras::k_sum(precision * (y_true - mean) ^ 2 + log_var, axis = 2)
  }


  spec_mse <- function(y_true, y_pred) {
    mean <- y_pred[, 1:output_dim]

    keras::k_log(keras::k_mean((y_true - mean) ^ 2))/2.303
  }

  model %>% keras::compile(
    optimizer = "adam",
    loss = heteroscedastic_loss,
    metrics = c(keras::custom_metric("spec_mse", spec_mse))
  )
  model
}







#' Create neural network with estimation of uncertainty through dropout
#'
#' @param combined_frame
#' @param output predicted value, one of "surv", "pseudo_eras", "expan"
#' @param transformation mathematical transformation applied to predicted value (function)
#'
#' @export
#'

create_b_model <- function(combined_frame, output = "surv", transformation = NULL, model = NULL) {
  indexes = sample(1:dim(combined_frame)[1], size = round(dim(combined_frame)[1]*0.85))
  predictors <- as.data.frame(combined_frame %>% dplyr::select(-surv, -pseudo_eras, -expan))
  predictors <- as.matrix(predictors)

  if (is.null(transformation)) {
    target <- as.matrix(combined_frame[output] , ncol = 1)
  }else{
    target <- as.matrix(transformation(combined_frame[output]) , ncol = 1)

  }


  xtrain = predictors[indexes,]
  indexes_n_val = sample(1:length(indexes), size = round(length(indexes)*0.85))
  xval = xtrain[-indexes_n_val,]
  xtrain = xtrain[indexes_n_val,]

  xtest = predictors[-indexes,]

  ytrain = target[indexes]
  yval = ytrain[-indexes_n_val]
  ytrain = ytrain[indexes_n_val]
  ytest = target[-indexes]


  # sample size (training data)
  n_train <- length(indexes_n_val)
  # prior length-scale
  input_dim <- dim(predictors)[2]
  output_dim <- 1
  # bayesian dropout model
  if (is.null(model)) {
    model <- model_b(n_train = n_train, input_dim = input_dim, output_dim = output_dim)
  }

  drop_track <- DropTracker$new()

  history <- model %>% keras::fit(x = xtrain, y = ytrain, epochs = 500, verbose = 0, validation_data = list(xval, yval),
                                  callbacks = list(
                                    keras::callback_early_stopping(patience = 70, restore_best_weights = T),
                                    keras::callback_reduce_lr_on_plateau(factor = 0.5, patience = 10),
                                    drop_track)
  )

  ps <- drop_track$ps %>% dplyr::as_tibble() %>% dplyr::select_if(function(x){sum(is.na(x)) == 0}) %>%
    dplyr::rename_with(~gsub("V", "layer", .x)) %>% tibble::rowid_to_column("epoch")
  print(ps)
  ps <- ps %>% tidyr::pivot_longer(cols = dplyr::starts_with("layer"))

  history$ps <- ps
  return( list(model = model, history = history, train_n_val_indexes = indexes, train_indexes = indexes_n_val))
}




DropTracker <- R6::R6Class("LossHistory",
                           inherit = keras::KerasCallback,

                           public = list(

                             ps  = NULL,

                             on_epoch_end = function(epoch, logs = list()) {
                               vector <- sapply(self$model$layers, FUN = function(x){
                                 v <- try(x$p_logit, silent = T)
                                 if (sum(class(v) == "try-error") == 1) {
                                   return(NA)
                                 }else{
                                   return(as.numeric(keras::k_sigmoid(v)))
                                 }
                               })

                               self$ps <- rbind(self$ps, vector)
                             }
                           ))








#' Evaluate model performance
#'
#' @param combined_frame
#' @param c_model_res
#' @param output
#' @param transformation
#' @export
evaluate_b_model <- function(combined_frame, c_model_res, output = "surv", transformation = NULL){
  # traing history (loss is heteroscedastic loss, (repetition))
  metrics <- dplyr::as_tibble(c_model_res$history$metrics)
  metrics$epoch <- 1:dim(metrics)[1]
  metrics <- metrics %>% tidyr::pivot_longer(cols = !epoch, names_to = "names")
  metrics <- metrics %>% dplyr::mutate(validation = grepl(names, pattern =  "val"),
                                       names = gsub(names, pattern =  "val_", replacement =  ""))

  history_plot <- metrics %>% ggplot2::ggplot(ggplot2::aes(epoch, value, col = validation)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::facet_grid(names~., scales = "free_y")

  dropout_freq_plot <- ggplot2::ggplot(data = c_model_res$history$ps, ggplot2::aes(epoch, value, color = name)) +
    ggplot2::geom_line(size = 1) + ggplot2::scale_y_log10()



  predictors <- as.data.frame(combined_frame %>% dplyr::select(-surv, -pseudo_eras, -expan))
  predictors <- as.matrix(predictors)


  num_MC_samples <- 100
  train_index <- 1:dim(combined_frame)[1] %in% c_model_res$train_n_val_indexes[c_model_res$train_indexes]

  MC_samples <- array(0, dim = c(num_MC_samples, length(train_index), 2 * 1))
  for (k in 1:num_MC_samples) {
    MC_samples[k, , ] <- (c_model_res$model %>% stats::predict(  as.matrix(combined_frame[,-(1:3)])  ))
  }
  # the means are in the first output column
  means <- MC_samples[, , 1]
  # average over the MC samples
  predictive_mean <- apply(means, 2, mean)
  epistemic_uncertainty <- apply(means, 2, var) # model uncertainty
  output_dim <- 1

  logvar <- MC_samples[, , (output_dim + 1):(output_dim * 2)]
  aleatoric_uncertainty <- exp(colMeans(logvar))  # implicit data uncertainty


  plot_d <- data.frame(
    predicted = predictive_mean,
    epistemic_uncertainty = sqrt(epistemic_uncertainty),
    aleatoric_uncertainty = sqrt(aleatoric_uncertainty),
    total_uncertainty = sqrt(epistemic_uncertainty + aleatoric_uncertainty),
    group = train_index,
    simulated = combined_frame[,output]
  )

  if (!is.null(transformation)) {
    plot_d <- plot_d %>% dplyr::mutate(simulated = transformation(simulated) )
  }

  metrics <- plot_d %>% dplyr::filter(!group) %>% dplyr::summarise(mse = mean( (predicted - simulated)^2 ),
                                                                   mae = mean( abs(predicted - simulated) ),
                                                                   rsq = (stats::cor(predicted, simulated))^2 )

  # regression plot
  reg_plot <- plot_d %>%
    ggplot2::ggplot(ggplot2::aes(simulated, predicted)) + ggplot2::geom_point(ggplot2::aes(col = group) ,alpha = 0.25, shape = 16 ) +
    ggplot2::geom_smooth(ggplot2::aes(group = as.factor(group))) + ggplot2::labs(color = "Used for\ntraining") +
    ggplot2::geom_abline(slope = 1, col = "red", size = 1) + ggplot2::theme_bw()

  # error histogram plot
  error_plot <- plot_d %>% ggplot2::ggplot(ggplot2::aes(  (predicted - simulated))) +
    ggplot2::geom_histogram(ggplot2::aes(fill = group, y = ggplot2::after_stat(ncount)), position = "identity", alpha = 0.5) +
    ggplot2::labs(fill = "Used for\ntraining") + ggplot2::theme_bw()

  # plot of total uncertainty ~ predicted value
  uncertainty_plot <- plot_d %>%
    ggplot2::ggplot(ggplot2::aes(simulated, total_uncertainty )) +
    ggplot2::geom_point(ggplot2::aes(col = group) ,alpha = 0.25, shape = 16 ) +
    ggplot2::geom_smooth(ggplot2::aes(group = as.factor(group))) + ggplot2::labs(color = "Used for\ntraining") +
    ggplot2::ylim(c(0,NA))

  return(list(metrics = metrics, predictions = plot_d,
              regression_plot = reg_plot, error_plot = error_plot, uncertainty_plot = uncertainty_plot,
              history_plot = history_plot, drop_plot = dropout_freq_plot))
}




#' Get prediction for new data
#'
#' @param assortativities_original data compiled using compile_assortativity_data function
#' @param simulation_data data used for model creation. Output of preprocess_simulation function
#' @param model_bayes generated model
#'
#' @return model predictions for new data
#' @export
get_predictions <- function(assortativities_original, simulation_data, model_bayes){
  if (!"matrix" %in% class(assortativities_original)) {
    assortativities_original <- as.matrix(assortativities_original)
  }

  s_par <- simulation_data$scale_params
  predictors <- (assortativities_original - s_par$`scaled:center`) / s_par$`scaled:scale`
  predictors <- predictors %*% simulation_data$rot_matrix

  prediction <- replicate(200, {model_bayes$model %>% predict(predictors)})

  predictive_mean <- mean(prediction[1,1,])
  epistemic_uncertainty <- var(prediction[1,1,])
  aleatoric_uncertainty <- exp(mean(prediction[1,2,]))

  prediction <- c(
    prediction = predictive_mean,
    epistemic_uncertainty = sqrt(epistemic_uncertainty),
    aleatoric_uncertainty = sqrt(aleatoric_uncertainty),
    u_overall = sqrt(epistemic_uncertainty) + sqrt(aleatoric_uncertainty)
  )
  prediction
}




preprocess_simulation <- function(simmulation_data, init_size){
  params_df <- purrr::map_dfr(simmulation_data, function(x){
    s <- x[[2]][1]
    e <- x[[2]][2]
    list(surv =  s/init_size, pseudo_eras = log2((init_size - s)/e + 1),  expan = e/s )
  })
  # for double assortativity
  assorts_data <- compile_assortativity_data(simmulation_data)
  assorts_data <- as.matrix(assorts_data)
  scaled <- scale(assorts_data)
  scale_params <- attributes(scaled)
  rot_mat <- prcomp(scaled, rank. = 20, scale. = T)$rot
  predictors <- scaled %*% rot_mat
  comb_frame <- cbind(params_df ,predictors)
  return(list(combined_frame = comb_frame, rot_matrix = rot_mat, scale_params = scale_params))
}




#' Predict regeneration for a given graph
#'
#' @param graph undirected graph with weight edge attribute and colour node attribute
#' @param survival_intervals range of values of survival fraction for simulation. If it is a single value, simulation and prediction of survival fraction is skipped and only expansion fraction is predicted
#' @param rel_expansion_intervals range of values of expansion fraction for simulation
#' @param with_push implement repulsing force in simulation?
#' @param alpha damping factor/factors
#'
#' @return returns list containing predictions, intermediate data, models and their evaluation. Prediction for expansion is square root of actual value due to improved prediction.
#' @export
get_complete_results_bayes <- function(graph, survival_intervals = c(0.5, 1), rel_expansion_intervals = c(0,1),
                                       with_push = T, alpha = c(0.4, 0.7),
                                       dsize_surv = 4000, dsize_expan = 3000){
  # values for original data
  ass_ori_data <- purrr::map(alpha, function(x){
    val <- assortativity_local_par(graph, igraph::get.vertex.attribute(graph)$colour, alpha = x)
    return(val)
  })
  assorts_data_ori <- compile_assortativity_data(list(list(ass_ori_data)))

  init_size <- igraph::vcount(graph)

  if (length(survival_intervals) != 1) {
    # simulation for survival
    set.seed(42)
    survivor_fraction <- stats::runif(dsize_surv, survival_intervals[1], survival_intervals[2])
    survivor_fraction[survivor_fraction > 1] <- 1
    set.seed(43)
    expanding_fraction <-  stats::runif(dsize_surv, rel_expansion_intervals[1], rel_expansion_intervals[2])
    expanding_fraction <- (1 - survivor_fraction)*1.1*expanding_fraction/survivor_fraction
    expanding_fraction[expanding_fraction > 1] <- 1
    sim_params <- cbind(survivor_fraction, expanding_fraction)
    set.seed(123)
    res_surv <- generate_data_fast(graph, sim_params, alpha = alpha, n_try = 1, with_push = with_push)
    surv_data <- preprocess_simulation(res_surv, init_size)


    # create , evaluate and use ML model
    model_surv_bayes <- reticulate::py_suppress_warnings( # each call raises multiple warnings about missing gradients for no reason
      create_b_model(surv_data$combined_frame)
      )
    eval_surv <- evaluate_b_model(surv_data$combined_frame, model_surv_bayes)
    survival_prediction <- get_predictions(assorts_data_ori, surv_data, model_surv_bayes)
    survival_fraction <- survival_prediction["prediction"]
    names(survival_prediction) <- paste0("surv_", names(survival_prediction))

    bhattach_dist <- get_Bhattacharyya_dist(ass_ori_data, res_surv, which_ass = 2)
    bhatt_tib <- surv_data$combined_frame[,1:3] %>% tibble::add_column(dist = bhattach_dist)
    bhatt_surv <- bhatt_tib %>% ggplot2::ggplot(ggplot2::aes(surv, dist)) +
      ggplot2::annotate("rect",xmin = survival_prediction[1] - survival_prediction[4],
                        xmax = survival_prediction[1] + survival_prediction[4],
                        ymin = -Inf, ymax = Inf, alpha = 0.25, fill = "red") +
      ggplot2::geom_point(size = 2, alpha = 0.2, shape = 16) +
      ggplot2::geom_vline(xintercept = survival_prediction[1], col = "red")
  } else{
    survival_fraction <- survival_intervals
  }



  # simulation for survival
  survivor_fraction <- round(survival_fraction * init_size, 0) / init_size
  survivor_fraction <- rep(survivor_fraction, dsize_expan)
  set.seed(42)
  expanding_fraction <-  runif(dsize_expan, 1/init_size, ((1 - survivor_fraction)*1.1/survivor_fraction)^0.5)^2
  sim_expansion <- cbind(survivor_fraction, expanding_fraction)
  set.seed(43)
  res_expansion <- generate_data_fast(graph, sim_expansion, alpha = alpha, n_try = 1, with_push = with_push)
  expan_data <- preprocess_simulation(res_expansion, init_size)


  # create , evaluate and use ML model
  model_expan_bayes <- reticulate::py_suppress_warnings( # each call raises multiple warnings about missing gradients
    create_b_model(expan_data$combined_frame, output = "expan", transformation = sqrt )
    )
  eval_expan <- evaluate_b_model(expan_data$combined_frame, model_expan_bayes, output = "expan", transformation = sqrt )
  expan_prediction <- get_predictions(assorts_data_ori, expan_data, model_expan_bayes)
  names(expan_prediction) <- paste0("expan_", names(expan_prediction))


  bhattach_dist <- get_Bhattacharyya_dist(ass_ori_data, res_expansion, which_ass = 2)
  bhatt_tib <- expan_data$combined_frame[,1:3] %>% tibble::add_column(dist = bhattach_dist)

  bhatt_expan <- bhatt_tib %>% ggplot2::ggplot(ggplot2::aes(expan, dist)) +
    ggplot2::annotate("rect",xmin = (expan_prediction[1] - expan_prediction[4])^2,
                      xmax = (expan_prediction[1] + expan_prediction[4])^2,
                      ymin = -Inf, ymax = Inf, alpha = 0.25, fill = "red") +
    ggplot2::geom_point(size = 2, alpha = 0.2, shape = 16) +
    ggplot2::geom_vline(xintercept = expan_prediction[1]^2, col = "red")


  if (length(survival_intervals) != 1) {
    return(list(survival = list(res = res_surv, data = surv_data, model_data = model_surv_bayes,
                                evaluation = eval_surv),
                expansion = list(res = res_expansion, data = expan_data, model_data = model_expan_bayes,
                                 evaluation = eval_expan),
                graph_data = list(assortativities = ass_ori_data, data = assorts_data_ori,
                                  # predictions
                                  predictions = c(survival_prediction, expan_prediction), size = init_size,
                                  survival_bhattacharyya = bhatt_surv, expansion_bhattacharyya = bhatt_expan)
    )
    )
  }else{
    return(list(expansion = list(res = res_expansion, data = expan_data, model_data = model_expan_bayes,
                                 evaluation = eval_expan),
                graph_data = list(assortativities = ass_ori_data, data = assorts_data_ori,
                                  # predictions
                                  predictions = expan_prediction, size = init_size,
                                  expansion_bhattacharyya = bhatt_expan)
    ))
  }
}





#' Small sample graph
#'
#' Small graph with 355 vertices, fragment of vessels of mouse 21 days after irradiation. Used for testing package installation
#'
#' @format Igraph graph. Important fr functionality are vertex attribute "colour", and edge attribute "weights"
"test_graph"





#' Test package functionality
#'
#' Run core function, get_complete_results on a test graph included in package. prediction quality is reduced, by limiting number simulations.
#'
#' @return result of get_complete_results
#'
#' @examples
#' @export
test_functionality <- function() {
  get_complete_results_bayes(RainbowGraph::test_graph, dsize_surv = 1000, dsize_expan = 1000)
}


