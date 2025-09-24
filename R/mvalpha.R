
#' Estimate Multi-Valued Krippendorff's Alpha
#'
#' `mvalpha()` calculates Krippendorff's alpha statistic when multi-valued observers are
#' allowed to apply multiple values to an observation.
#'
#' @param data a data frame containing a list column for each observer. Each row represents
#' an observation unit, and each cell contains a vector of 0 to w unique values, where w is the
#' number of unique values found in the data set. `NA` values are used to represent
#' missing observations and `NULL` values represent the empty set, `{}`, of responses.
#' @param type a string describing the data type of the label set. This can be "nominal",
#' "ordinal", "interval", or "ratio" and is used to select the appropriate distance metric.
#' @param verbose a logical value which toggles whether status updates are printed to
#' the console while alpha is being calculated.
#' @param n_boot an integer representing the number of bootstrap estimates to calculate
#' for mvDo. The default, `NULL`, will not generate additional estimates.
#' @param parallelize a logical value indicating whether to implement parallelization
#' using the `parallel` package.
#' @param cluster_size an integer describing the number of cores to allocate to parallelization.
#' If `NULL` and `parallelize=TRUE`, then the maximum number of available cores minus 1 will
#' be used.
#' @returns An object of class `mvalpha`
#' @export
#' @example man/examples/mvalpha_example.R
#' @references
#' \insertRef{Krippendorff-Craggs-2016}{mvalpha}
#' @import parallel
#' @importFrom Rdpack reprompt
#' @importFrom stats na.omit setNames
#' @importFrom utils combn head tail


mvalpha <-
  function(data, type = "nominal", verbose = TRUE, n_boot = NULL, parallelize = FALSE, cluster_size = NULL){

    calc_prods_and_metrics <-
      function(C, `K_|K|`){
        lapply(`K_|K|`, function(K){
          inner_set_ops <- set_ops(K, C, type = type)
          prod_all <- prod(n_c_w[c(C, inner_set_ops$A_diff_B) + null_set_observed], na.rm = TRUE)
          prod_int <- prod(n_c_w[inner_set_ops$A_intersect_B + null_set_observed] - 1, na.rm = TRUE)
          metric <- metric_delta_CK(C, K, inner_set_ops)
          c(prod_all * prod_int, metric)
        }) |> do.call(rbind, args = _) # calculate over K of cardinality |K|
      }

    bootstrap_mvDo <-
      function(){
        vapply(units[which(m != 0)], function(u){
          (sample(dist_CK, size = (m[u] - 1) * m[u] / 2, replace = TRUE, prob = p_CK) /
             (m[u] - 1)) |>
            sum()
        }, 1) |> sum() * 2 / n..
      }

    metric_delta_CK <-
      function(C, K, KC_set_ops){
        if(type == "nominal"){
          if(!length(C) & !length(K)){return(0)}
          else{
            result <-
              1 - (2 * length(KC_set_ops$A_intersect_B) /
                     (length(C) + length(K)))
            return(result)
          }
        }
        if(rlang::is_empty(C) | rlang::is_empty(KC_set_ops$A_diff_B)){lhs_numerator <- 0}else{
          lhs_numerator <-
            outer(C, KC_set_ops$A_diff_B, metric_delsq_ck) |>
            unlist() |>
            sum(na.rm = TRUE) |>
            (\(x) x / length(C))()
        }
        if(rlang::is_empty(K) | rlang::is_empty(KC_set_ops$B_diff_A)){rhs_numerator <- 0}else{
          rhs_numerator <-
            outer(K, KC_set_ops$B_diff_A, metric_delsq_ck) |>
            unlist() |>
            sum(na.rm = TRUE) |>
            (\(x) x / length(K))()
        }
        denominator <- length(C) + length(K)
        result <-
          ifelse(lhs_numerator == 0 & rhs_numerator == 0 & denominator == 0,
                 0,
                 ifelse(rlang::is_empty(C) | rlang::is_empty(K),
                        1,
                        (lhs_numerator + rhs_numerator) / denominator))
        return(result)
      }

    nominal_delsq_ck <-
      function(c, k){as.numeric(!(c == k))}

    ordinal_delsq_ck <-
      Vectorize(
        function(c, k){
          if(is.character(c)) c <- match(values[c], labels)
          if(is.character(k)) k <- match(values[k], labels)
          ((sum(n_c_w[c:k]) - ((n_c_w[c] + n_c_w[k]) / 2)) /
              (n.. - (n_c_last + n_c_first) / 2)) ^ 2
        },
        vectorize.args = c("c", "k")
      )

    interval_delsq_ck <-
      function(c, k){
        ((c - k) / (c_max - c_min)) ^ 2
      }

    ratio_delsq_ck <-
      function(c, k){
        ((c - k) / (c + k)) ^ 2
      }

    # Select appropriate metric

    metric_delsq_ck <-
      switch(type,
             nominal = nominal_delsq_ck,
             ordinal = ordinal_delsq_ck,
             interval = interval_delsq_ck,
             ratio = ratio_delsq_ck)

    # Set up cluster if parallelization is enabled

    if(parallelize){
      if(verbose){cat("Initializing cluster\n")}
      if(is.null(cluster_size)){cluster_size <- parallel::detectCores() - 1}
      cl <- parallel::makeCluster(cluster_size)
    }

    # Summary vectors and tables used to organize data

    labels <- data |> unlist() |> unique() |> sort()
    observers <- colnames(data)
    units <- rownames(data)
    n_labels <- length(labels)
    n_observers <- length(observers)
    n_units <- length(units)
    continuous_data <- type %in% c("interval", "ratio") # logical

    label_array <- # indicator array for labels by observer and unit
      apply(data, MARGIN = c(1, 2),
            function(x){
              if(!is.na(x)){as.numeric(labels %in% unlist(x))
              }else{
                rep(NA, length(labels))
              }
            }) |>
      aperm(c(2, 3, 1))

    dimnames(label_array) <- list(unit = units, observer = observers, label = labels)

    values_def <- apply(label_array, 3, rbind) |> unique(MARGIN = 1) |> stats::na.omit()
    n_values <- nrow(values_def)

    values <- apply(values_def, MARGIN = 1, function(x){
      if(type == "nominal"){x[which(x == 1)] |> names()}
      else{x[which(x == 1)] |> names() |> as.numeric()}})

    value_names <- values |> lapply(function(x){paste0("{", paste0(x, collapse = ", "), "}")}) |> unlist() |> unname()
    names(values) <- value_names

    rownames(values_def) <- value_names

    values_by_unit <- # create value-by-unit matrix
      apply(values_def, MARGIN = 1, function(v){
        apply(label_array, MARGIN = 1, function(x){
          apply(x, MARGIN = 1, function(r){
            all(v==r)
          }) |> sum(na.rm = TRUE)
        })
      }) |> t() |> unname()

    dimnames(values_by_unit) <- list(value = value_names, unit = units)

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

    if(type == "ordinal"){ # define additional variables needed for ordinal metric
      n_c_first <- utils::head(n_c_w, n = 1)
      n_c_last <- utils::tail(n_c_w, n = 1)
    }

    if(type == "interval"){ # define additional variables needed for interval metric
      c_min <- min(labels)
      c_max <- max(labels)
    }

    P <- # probability of observing a label set of cardinality |C|
      vapply(unique_cardinalities, function(x){
        sum(n[which(value_cardinalities == x)]) / n..
      }, 1) |> stats::setNames(unique_cardinalities)

    all_label_combinations <- # generate all ways to choose |C| from n_labels
      lapply(unique_cardinalities,
        function(f){if(f == 0){utils::combn(1, 1)}else{utils::combn(n_labels, f)}}) |>
      stats::setNames(unique_cardinalities)

    n_label_combinations <- lapply(all_label_combinations, ncol) |> unlist()

    if(continuous_data){lapply(all_label_combinations, function(`|H|`){ # redefine indicator combinations using actual values if data are continuous
      matrix(labels[`|H|`], nrow = nrow(`|H|`))
    })}

    p_CK <- # probabilities of the observed coincidence of C-K pairs
      lapply(value_names, function(C){
        lapply(value_names, function(K){
          if(C == K){(values_by_unit[C, ] * (values_by_unit[K, ] - 1) / (m - 1)) |> sum()}
          else{(values_by_unit[C, ] * values_by_unit[K, ] / (m - 1)) |> sum()}
        }) |> unlist()
      }) |> do.call(cbind, args = _) |>
      (\(i) i / n..)()

    dimnames(p_CK) <- list(value_names, value_names)

    dist_CK <- # distances between C-K pairs
      lapply(values, function(C){
        lapply(values, function(K){
          metric_delta_CK(C, K, set_ops(K, C, type = type))
        }) |> unlist()
      }) |> do.call(cbind, args = _)

    # Print status update

    if(verbose){cat(
      paste0("n Units: ", n_units, "\n",
             "n Observers: ", n_observers, "\n",
             "n Labels: ", n_labels, "\n",
             "n Values: ", n_values, "\n",
             "Max Cardinality: ", max(unique_cardinalities), "\n",
             "Possible C-K pairs: ",
             (n_label_combinations %*% t(n_label_combinations)) |>
               (\(m) m[upper.tri(m, diag = TRUE)])() |>
               sum() |> prettyNum(big.mark = ","),
             "\n"
             ))}

    # Notation:
    # C is a set (vector) of elements
    # |C| is the cardinality (integer) of set C
    # C_|C| is the collection (list) of sets with cardinality |C|

    mvDo <- sum(p_CK * dist_CK)

    mvDe <-
      vapply(unique_cardinalities, function(`|C|`){
        `|C|_str` <- as.character(`|C|`)
        `C_|C|` <- if(`|C|` == 0){0}
        else{split(all_label_combinations[[`|C|_str`]], col(all_label_combinations[[`|C|_str`]]))}

        vapply(unique_cardinalities[which(unique_cardinalities >= `|C|`)], function(`|K|`){
          `|K|_str` <- as.character(`|K|`)
          `K_|K|` <- if(`|K|` == 0){0}
          else{split(all_label_combinations[[`|K|_str`]], col(all_label_combinations[[`|K|_str`]]))}

          if(verbose){cat(paste0( # Status update
            "|C| = ", `|C|`, ", |K| = ", `|K|`,
            ", C-K pairs to evaluate = ",
            prettyNum(n_label_combinations[as.character(`|C|`)] *
            n_label_combinations[as.character(`|K|`)], big.mark = ","),
            "\n"))}

          prods_and_metrics <- # Pi products and metrics
            if(parallelize){
              parallel::parLapply(cl, `C_|C|`, function(C){calc_prods_and_metrics(C, `K_|K|`)}) |>
                do.call(rbind, args = _)}
            else{lapply(`C_|C|`, function(C){
              calc_prods_and_metrics(C, `K_|K|`)})} |>
                do.call(rbind, args = _)

          P[`|C|_str`] * # Multiply by proportion of observed sets with cardinality |C|
            P[`|K|_str`] *  # Multiply by proportion of observed sets with cardinality |K|
            (sum(prods_and_metrics[, 1] * prods_and_metrics[, 2]) /
               sum(prods_and_metrics[, 1]) *
                (1 + (`|C|_str` != `|K|_str`))) # Double when off main diagonal

        }, 1) |> sum(na.rm = TRUE) # sum over |K|
      }, 1) |> sum(na.rm = TRUE) # sum over |C|

    mvDo_boot <- # Bootstrap method described by Krippendorff and Craggs (2016)
      if(!is.null(n_boot)){
        if(parallelize){
          parallel::parLapply(cl, 1:n_boot, function(iter){bootstrap_mvDo()}) |>
            unlist()}
        else{lapply(1:n_boot, function(iter){bootstrap_mvDo()}) |>
            unlist()}}
      else{NULL}

    if(parallelize){parallel::stopCluster(cl)}

    return(
      new_mvalpha(
        mvalpha = 1 - mvDo / mvDe,
        type = type,
        mvDo = mvDo,
        mvDe = mvDe,
        bootstrap_mvalpha = 1 - (mvDo_boot / mvDe),
        unique_cardinalities = unique_cardinalities,
        units = units,
        observers = observers,
        labels = labels,
        values = values,
        values_by_unit = values_by_unit,
        dist_CK = dist_CK,
        p_CK = p_CK,
        data = data
      )
    )

  }







