
#' @export
#' @import parallel


mv_alpha <-
  function(data, type = "nominal", verbose = TRUE, n_boot = NULL, parallelize = FALSE, cluster_size = NULL){

    # Internal functions

    calc_prods_and_metrics <-
      function(C, `K_|K|`){
        C <- unlist(C)
        lapply(`K_|K|`, function(K){
          K <- unlist(K)
          inner_set_ops <- mvalpha::set_ops(K, C)
          prod_ns <-
            prod(n_c_w[C + null_set_observed], na.rm = TRUE) * # prod c in C
            prod(n_c_w[inner_set_ops$A_diff_B + null_set_observed], na.rm = TRUE) * # prod unique k in K
            prod(n_c_w[inner_set_ops$A_intersect_B + null_set_observed] - 1, na.rm = TRUE) # prod intersect C and K
          metric <- mvalpha::metric_delta_CK(C, K, inner_set_ops, type = type)
          c(prod_ns, metric)
        }) |> do.call(rbind, args = _) # calculate over K of cardinality |K|
      }

    bootstrap_mvDo <-
      function(){
        vapply(1:n_units, function(u){
          (sample(dist_CK, size = (m[u] - 1) * m[u] / 2, replace = TRUE, prob = p_CK) /
             (m[u] - 1)) |>
            sum()
        }, 1) |> sum() * 2 / n..
      }

    if(parallelize){
      if(verbose){cat("Initializing cluster\n")}
        if(is.null(cluster_size)){cluster_size <- parallel::detectCores() - 1}
        cl <- parallel::makeCluster(cluster_size)
      }

    features <- data |> unlist() |> unique() |> sort()
    readers <- colnames(data)
    units <- rownames(data)
    n_features <- length(features)
    n_readers <- length(readers)
    n_units <- length(units)

    feature_array <- # indicator array for features by reader and unit
      apply(data, MARGIN = c(1, 2),
            function(x){
              if(!is.na(x)){as.numeric(features %in% unlist(x))
              }else{
                rep(NA, length(features))
              }
            }) |>
      aperm(c(2, 3, 1))

    dimnames(feature_array) <- list(unit = units, reader = readers, feature = features)

    values_def <- apply(feature_array, 3, rbind) |> unique(MARGIN = 1) |> na.omit()
    n_values <- nrow(values_def)

    values <- apply(values_def, MARGIN = 1, function(x){x[which(x == 1)] |> names()})
    value_names <- values |> lapply(function(x){paste0("{", paste0(x, collapse = ", "), "}")}) |> unlist() |> unname()
    names(values) <- value_names

    rownames(values_def) <- value_names

    values_by_unit <- # create value by unit matrix
      apply(values_def, MARGIN = 1, function(v){
        apply(feature_array, MARGIN = 1, function(x){
          apply(x, MARGIN = 1, function(r){
            all(v==r)
          }) |> sum(na.rm = TRUE)
        })
      }) |> t()

    m <- colSums(values_by_unit)
    m[which(m<2)] <- 0 # remove non-pairable observations
    values_by_unit[, which(m==0)] <- 0

    n <- rowSums(values_by_unit)
    n.. <- sum(n)

    value_cardinalities <- rowSums(values_def)
    unique_cardinalities <- unique(value_cardinalities) |> sort()
    null_set_observed <- 0 %in% unique_cardinalities # special behavior flag for null set

    w <- data |> unlist() |> unique() |> setdiff(NA) |> length()
    if(null_set_observed) w <- w + 1
    n_c_w <- colSums(values_def * n)
    if(null_set_observed) n_c_w <- c(null_set = 1, n_c_w)

    P <- # prob of observing a label set of cardinality |C|
      vapply(unique_cardinalities, function(x){
        sum(n[which(value_cardinalities == x)]) / n..
      }, 1) |> stats::setNames(unique_cardinalities)

    all_feature_combinations <- # all ways to choose |C| from n_features
      lapply(unique_cardinalities,
        function(f){
          if(f == 0){"null_set"
          }else{combn(n_features, f)}
        }) |> stats::setNames(unique_cardinalities)

    n_feature_combinations <- lapply(all_feature_combinations, ncol) |> unlist()
    if(null_set_observed) n_feature_combinations <- c(null_set = 1, n_feature_combinations)

    p_CK <- # probabilities of the observed coincidence of C-K pairs
      lapply(value_names, function(C){
        lapply(value_names, function(K){
          if(C == K){
            (values_by_unit[C, ] * (values_by_unit[K, ] - 1) / (m - 1)) |> sum()
          }else{
              (values_by_unit[C, ] * values_by_unit[K, ] / (m - 1)) |> sum()
            }
        }) |> unlist()
      }) |> do.call(cbind, args = _) |>
      (\(i) i / n..)()

    dimnames(p_CK) <- list(value_names, value_names)

    dist_CK <- # distances between C-K pairs
      lapply(values, function(C){
        lapply(values, function(K){
          mvalpha::metric_delta_CK(C, K, mvalpha::set_ops(K, C), type = type)
        }) |> unlist()
      }) |> do.call(cbind, args = _)

    if(verbose){cat(
      paste0("n Units: ", n_units, "\n",
             "n Readers: ", n_readers, "\n",
             "n Features: ", n_features, "\n",
             "n Values: ", n_values, "\n",
             "Max Cardinality: ", max(unique_cardinalities), "\n",
             "Possible C-K pairs: ", sum(n_feature_combinations %*% t(n_feature_combinations)), "\n"
             ))}

    # C is a set (vector) of elements
    # K is a set (vector) of elements
    # |C| is the cardinality (integer) of set C
    # |K| is the cardinality (integer) of set K
    # C_|C| is a list of sets (vectors) with cardinality |C|
    # K_|K| is a list of sets (vectors) with cardinality |K|

    mvDo <- sum(p_CK * dist_CK)

    mvDe <-
      vapply(unique_cardinalities, function(`|C|`){
        `C_|C|_str` <- as.character(`|C|`)
        `C_|C|` <- if(`|C|` == 0)
          {0}else{
            split(all_feature_combinations[[`C_|C|_str`]], col(all_feature_combinations[[`C_|C|_str`]]))}

        P[`C_|C|_str`] * # Multiply by proportion of observed sets with cardinality |C|
          vapply(unique_cardinalities, function(`|K|`){
            `K_|K|_str` <- as.character(`|K|`)
            `K_|K|` <- if(`|K|` == 0)
            {0}else{
              split(all_feature_combinations[[`K_|K|_str`]], col(all_feature_combinations[[`K_|K|_str`]]))}

            if(verbose){
              cat(
                paste0("|C| = ", `|C|`, ", |K| = ", `|K|`,
                       ", C-K pairs to evaluate = ",
                       n_feature_combinations[`|C|` + null_set_observed] * n_feature_combinations[`|K|` + null_set_observed], "\n"))
            }

            prods_and_metrics <- # Pi products and metrics
              if(parallelize){
                parallel::parLapply(cl, `C_|C|`, function(C){
                  calc_prods_and_metrics(C, `K_|K|`)
                  }) |> do.call(rbind, args = _)
              }else{
                lapply(`C_|C|`, function(C){
                  calc_prods_and_metrics(C, `K_|K|`)})
              } |> do.call(rbind, args = _)

            P[`K_|K|_str`] *  # Multiply by proportion of observed sets with cardinality |K|
              sum(prods_and_metrics[, 1] * prods_and_metrics[, 2]) /
              sum(prods_and_metrics[, 1])

          }, 1) |> sum(na.rm = TRUE) # sum over |K|
      }, 1) |> sum(na.rm = TRUE) # sum over |C|

    mvDo_boot <-
      if(!is.null(n_boot)){
        if(parallelize){
          parallel::parLapply(cl, 1:n_boot, function(iter){
            bootstrap_mvDo()
          }) |> unlist()
        }else{
          lapply(1:n_boot, function(iter){
            bootstrap_mvDo()
          }) |> unlist()
        }
      }else{
        NULL
      }

    if(parallelize){parallel::stopCluster(cl)}

    return(list(
      mvalpha = 1 - mvDo / mvDe,
      bootstrap_mvalpha = 1 - (mvDo_boot / mvDe),
      mvDo = mvDo,
      mvDe = mvDe,
      unique_cardinalities = unique_cardinalities,
      features = features,
      units = units,
      readers = readers,
      data = data
    ))

  }







