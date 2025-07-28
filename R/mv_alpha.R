
#' @export
#' @import parallel

mv_alpha <-
  function(data, type = "nominal", clusters = parallel::detectCores() - 1){

    cl <- parallel::makeCluster(clusters)

    parallel::clusterExport(cl, list("set_ops", "metric_delsq_ck", "metric_delta_CK"))

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
                            }else{
                              combn(n_features, f)
                            }
                          }
      ) |> stats::setNames(unique_cardinalities)

    # C is a set (vector) of elements
    # K is a set (vector) of elements
    # |C| is an integer
    # |K| is an integer
    # C_|C| is a list of sets (vectors) of elements
    # K_|K| is a list of sets (vectors) of elements

    mvDo <-
      vapply(units[which(m > 0)], function(u){
        unit_values <- value_names[which(values_by_unit[, u] > 0)]
        lapply(unit_values, function(c_ind){
          values_by_unit[c_ind, u] *
            lapply(unit_values, function(k){
              values_by_unit[k, u] * metric_delta_CK(values[[c_ind]], values[[k]])
            }) |> unlist() |> sum(na.rm = TRUE) # sum over k
        }) |> unlist() |> sum(na.rm = TRUE) |> # sum over c
          (\(x) x / (m[u] - 1))()
      }, 1) |> sum(na.rm = TRUE) |> # sum over u
      (\(x) x / n..)()

    mvDe <-
      vapply(unique_cardinalities, function(`|C|`){
        `C_|C|_str` <- as.character(`|C|`)
        `C_|C|` <- if(`|C|` == 0)
          {0}else{
          split(all_feature_combinations[[`C_|C|_str`]], col(all_feature_combinations[[`C_|C|_str`]]))}
        print(`|C|`)
        P[`C_|C|_str`] * # Multiply by proportion of observed sets with cardinality |C|
          vapply(unique_cardinalities, function(`|K|`){
            `K_|K|_str` <- as.character(`|K|`)
            `K_|K|` <- if(`|K|` == 0)
            {0}else{
              split(all_feature_combinations[[`K_|K|_str`]], col(all_feature_combinations[[`K_|K|_str`]]))}
            print(`|K|`)
            denominator <-
              parallel::parLapply(cl, `C_|C|`, function(C_inner)
              {
                C_inner <- unlist(C_inner)
                lapply(`K_|K|`, function(K_inner)
                {
                  K_inner <- unlist(K_inner)

                  inner_set_ops <- set_ops(K_inner, C_inner)

                  prod(n_c_w[C_inner + null_set_observed], na.rm = TRUE) * # prod c in C
                    prod(n_c_w[inner_set_ops$A_diff_B + null_set_observed], na.rm = TRUE) * # prod unique k in K
                    prod(n_c_w[inner_set_ops$A_intersect_B + null_set_observed] - 1, na.rm = TRUE) # prod intersect C and K
                }) |> unlist() |> sum(na.rm = TRUE) # sum over K of cardinality |K|
              }) |> unlist() |> sum(na.rm = TRUE) # sum over C of cardinality |C|

            P[`K_|K|_str`] *  # Multiply by proportion of observed sets with cardinality |K|
              parallel::parLapply(cl, `C_|C|`, function(C_outer){
                C_outer <- unlist(C_outer)
                lapply(`K_|K|`, function(K_outer){
                  K_outer <- unlist(K_outer)

                  outer_set_ops <- set_ops(K_outer, C_outer)

                  numerator <-
                    prod(n_c_w[C_outer + null_set_observed], na.rm = TRUE) * # prod c in C
                    prod(n_c_w[outer_set_ops$A_diff_B + null_set_observed], na.rm = TRUE) * # prod unique k in K
                    prod(n_c_w[outer_set_ops$A_intersect_B + null_set_observed] - 1, na.rm = TRUE) # prod intersect C and K

                  numerator / denominator * metric_delta_CK(C_outer, K_outer, "nominal")

                }) |> unlist() |> sum(na.rm = TRUE) # sum over K_outer of cardinality |K|
              }) |> unlist() |> sum(na.rm = TRUE) # sum over C_outer of cardinality |C|
          }, 1) |> sum(na.rm = TRUE) # sum over |K|
      }, 1) |> sum(na.rm = TRUE) # sum over |C|

    result <- 1 - mvDo / mvDe

    parallel::stopCluster(cl)

    return(list(
      result = result,
      unique_cardinalities = unique_cardinalities,
      features = features,
      units = units,
      readers = readers
    ))

  }
