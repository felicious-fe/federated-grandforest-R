##' @export
grandforest_sum_models <- function(model1, model2, keep.inbag = FALSE, probability = FALSE, seed = NULL, num.threads = NULL, verbose=TRUE) {
  ## Check forest argument
  if (class(model1$forest) != "grandforest.forest") {
    stop("Error: Invalid class of input object.")
  } else {
    forest <- model1$forest
  }
  if (is.null(forest$dependent.varID) || is.null(forest$num.trees) ||
    is.null(forest$child.nodeIDs) || is.null(forest$split.varIDs) ||
    is.null(forest$split.values) || is.null(forest$independent.variable.names) ||
    is.null(forest$treetype)) {
    stop("Error: Invalid forest object.")
  }
  if (forest$treetype == "Survival" && (is.null(forest$status.varID)  ||
    is.null(forest$chf) || is.null(forest$unique.death.times))) {
    stop("Error: Invalid forest object.")
  }

  ## Check forest argument
  if (class(model2$forest) != "grandforest.forest") {
    stop("Error: Invalid class of input object.")
  } else {
    forest <- model2$forest
  }
  if (is.null(forest$dependent.varID) || is.null(forest$num.trees) ||
    is.null(forest$child.nodeIDs) || is.null(forest$split.varIDs) ||
    is.null(forest$split.values) || is.null(forest$independent.variable.names) ||
    is.null(forest$treetype)) {
    stop("Error: Invalid forest object.")
  }
  if (forest$treetype == "Survival" && (is.null(forest$status.varID)  ||
    is.null(forest$chf) || is.null(forest$unique.death.times))) {
    stop("Error: Invalid forest object.")
  }


  if (model1$treetype == "Classification") {
    treetype <- 1
  } else if (model1$treetype == "Regression") {
    treetype <- 3
  } else if (model1$treetype == "Survival") {
    treetype <- 5
  } else if (model1$treetype == "Probability estimation") {
    treetype <- 9
  }
  print(paste(model1$treetype, "->", treetype))

  ## Keep inbag
  if (!is.logical(keep.inbag)) {
    stop("Error: Invalid value for keep.inbag")
  }

  ## Importance mode
  if (is.null(model1$importance.mode) || model1$importance.mode == "none") {
    importance.mode <- 0
  } else if (model1$importance.mode == "impurity") {
    importance.mode <- 1
  } else if (model1$importance.mode == "impurity_corrected" || model1$importance.mode == "impurity_unbiased") {
    importance.mode <- 5
    if (!is.null(split.select.weights)) {
      stop("Corrected impurity importance not supported in combination with split.select.weights.")
    }
  } else if (model1$importance.mode == "permutation") {
    if (scale.permutation.importance) {
      importance.mode <- 2
    } else {
      importance.mode <- 3
    }
  } else {
    stop("Error: Unknown importance mode.")
  }

  ## Splitting rule
  if (is.null(model1$splitrule)) {
    if (treetype == 5) {
      model1$splitrule <- "logrank"
    } else if (treetype == 3) {
      model1$splitrule <- "variance"
    } else if (treetype %in% c(1, 9)) {
      model1$splitrule <- "gini"
    }
    splitrule.num <- 1
  } else if (model1$splitrule == "logrank") {
    if (treetype == 5) {
      splitrule.num <- 1
    } else {
      stop("Error: logrank splitrule applicable to survival data only.")
    }
  } else if (model1$splitrule == "gini") {
    if (treetype %in% c(1, 9)) {
      splitrule.num <- 1
    } else {
      stop("Error: Gini splitrule applicable to classification data only.")
    }
  } else if (model1$splitrule == "variance") {
    if (treetype == 3) {
      splitrule.num <- 1
    } else {
      stop("Error: variance splitrule applicable to regression data only.")
    }
  } else if (model1$splitrule == "auc" || model1$splitrule == "C") {
    if (treetype == 5) {
      splitrule.num <- 2
    } else {
      stop("Error: C index splitrule applicable to survival data only.")
    }
  } else if (model1$splitrule == "auc_ignore_ties" || model1$splitrule == "C_ignore_ties") {
    if (treetype == 5) {
      splitrule.num <- 3
    } else {
      stop("Error: C index splitrule applicable to survival data only.")
    }
  } else if (model1$splitrule == "maxstat") {
    if (treetype == 5 || treetype == 3) {
      splitrule.num <- 4
    } else {
      stop("Error: maxstat splitrule applicable to regression or survival data only.")
    }
  } else if (model1$splitrule == "extratrees") {
    splitrule.num <- 5
  } else {
    stop("Error: Unknown splitrule.")
  }

  ## Num threads
  ## Default 0 -> detect from system in C++.
  if (is.null(num.threads)) {
    num.threads <- 0
  } else if (!is.numeric(num.threads) || num.threads < 0) {
    stop("Error: Invalid value for num.threads")
  }

  ## Seed
  if (is.null(seed)) {
    seed <- runif(1 , 0, .Machine$integer.max)
  }

  # create dummy data matrices with double 0,0,0,... vectors as values
  data1 <- matrix(rep(0, len=(model1$num.samples * (length(model1$forest$independent.variable.names) + 1))), nrow=model1$num.samples)
  data2 <- matrix(rep(0, len=(model2$num.samples * (length(model2$forest$independent.variable.names) + 1))), nrow=model2$num.samples)
  variable.names1 <- model1$forest$independent.variable.names
  variable.names2 <- model2$forest$independent.variable.names
  combined.independent.variable.names <- c(model1$forest$independent.variable.names,
                                           setdiff(model2$forest$independent.variable.names,
                                                   model1$forest$independent.variable.names))

  num.trees1 <- model1$forest$num.trees
  num.trees2 <- model2$forest$num.trees

  ## Defaults for variables not needed
  graph <- matrix()
  dependent.variable.name <- "none"
  mtry <- 0
  subgraph <- 0
  min.node.size <- 0
  split.select.weights <- list(c(0, 0))
  use.split.select.weights <- FALSE
  always.split.variables <- c("0", "0")
  use.always.split.variables <- FALSE
  status.variable.name <- "status"
  prediction.mode <- TRUE
  sample.with.replacement <- TRUE
  probability <- FALSE
  unordered.factor.variables <- c("0", "0")
  use.unordered.factor.variables <- FALSE
  save.memory <- FALSE
  alpha <- 0
  minprop <- 0
  case.weights <- c(0, 0)
  use.case.weights <- FALSE
  predict.all <- FALSE
  keep.inbag <- FALSE
  sample.fraction <- 1
  holdout <- FALSE
  prediction.type <- 1
  num.random.splits <- 1
  random.root <- FALSE

  combined_model <- grandforest_sum_modelsCpp(model1, model2, treetype, dependent.variable.name, data1, data2,
                                              variable.names1, variable.names2, graph, mtry,
                                              num.trees1, num.trees2, verbose, seed, num.threads, importance.mode, subgraph,
                                              min.node.size, split.select.weights, use.split.select.weights,
                                              always.split.variables, use.always.split.variables, status.variable.name,
                                              prediction.mode, sample.with.replacement, probability, unordered.factor.variables,
                                              use.unordered.factor.variables, save.memory, splitrule.num, case.weights,
                                              use.case.weights, predict.all, keep.inbag, sample.fraction, alpha,
                                              minprop, holdout, prediction.type, num.random.splits, random.root)

  if (length(combined_model) == 0) {
    stop("User interrupt or internal error.")
  }

  ## Prepare results
  if (importance.mode != 0) {
    combined.independent.variable.names <- unique(c(model1$forest$independent.variable.names, model2$forest$independent.variable.names))
    combined.variable.importance <- setNames(
                                      mapply(
                                        function(x1,x2,x1samples,x2samples) {
                                          if(length(x1)==0){x1=0}
                                          if(length(x2)==0){x2=0}
                                          return(x1*x1samples + x2*x2samples)
                                        },
                                        model1$variable.importance[combined.independent.variable.names],
                                        model2$variable.importance[combined.independent.variable.names],
                                        x1samples=model1$num.samples,
                                        x2samples=model2$num.samples
                                      ),
                                      combined.independent.variable.names
                                    )
    relative.combined.variable.importance <- as.list(as.numeric(combined.variable.importance)/sum(as.numeric(combined.variable.importance)))
    names(relative.combined.variable.importance) <- names(combined.variable.importance)
    combined_model$variable.importance <- relative.combined.variable.importance
  }
  if(length(combined_model$variable.frequency) > 0) {
    #names(combined_model$variable.frequency) <- model1$forest$independent.variable.names
  }

  ## Splitrule
  combined_model$splitrule <- model1$splitrule

  ## Set treetype
  combined_model$treetype <- model1$treetype
  combined_model$min.node.size <- min(model1$min.node.size, model2$min.node.size)
  combined_model$mtry <- max(model1$mtry, model2$mtry)
  if (treetype == 3) {
    #TODO is this correct
    combined_model$r.squared <- mean(model1$r.squared, model2$r.squared)
  }
  combined_model$call <- sys.call()
  combined_model$importance.mode <- model1$importance.mode
  combined_model$num.samples <- model1$num.samples + model2$num.samples
  combined_model$replace <- (model1$replace & model2$replace)

  ## Write forest object
  if (length(model1$forest$levels) > 0) {
    combined_model$forest$levels <- c(model1$forest$levels, model2$forest$levels)
  }
  # TODO This could be risky if the two datasets don't have the same variables
  combined_model$forest$independent.variable.names <- model1$forest$independent.variable.names
  combined_model$forest$treetype <- combined_model$treetype
  class(combined_model$forest) <- "grandforest.forest"

  class(combined_model) <- "grandforest"
  return(combined_model)
}


symmetric_difference <- function(set1, set2) {
  symmetric_difference_vector <- c(setdiff(set1, set2), setdiff(set2, set1))
  return(unique(symmetric_difference_vector))
}