#' @export

set_ops <-
  function (A, B, type = "nominal") {

    # https://stackoverflow.com/a/72631719

    ind <- match(B, A, nomatch = 0)
    A_intersect_B <- A[ind]
    A_diff_B <- A[-c(ind, length(A) + 1)]  # https://stackoverflow.com/a/52772380

    if(type == "nominal"){B_diff_A <- NULL}
    else{
      ind2 <- match(A, B, nomatch = 0)
      B_diff_A <- B[-c(ind2, length(B) + 1)]
      }

    return(
      list(
        A_intersect_B = A_intersect_B,
        A_diff_B = A_diff_B,
        B_diff_A = B_diff_A
      )
    )
  }

#' @export

metric_delta_CK <-
  function(C, K, KC_set_ops, type = "nominal"){

    if(type == "nominal"){
      if(!length(C) & !length(K)){return(0)}
      else{
        result <-
          1 - (2 * length(KC_set_ops$A_intersect_B) /
            (length(C) + length(K)))
        return(result)
      }
    }

    if(rlang::is_empty(C)){lhs_numerator <- 0}else{
      lhs_numerator <-
        outer(C, KC_set_ops$A_diff_B, mvalpha::metric_delsq_ck, type = type) |>
          sum(na.rm = TRUE) |>
          (\(x) x / length(C))()
    }

    if(rlang::is_empty(K)){rhs_numerator <- 0}else{
      rhs_numerator <-
        outer(K, KC_set_ops$B_diff_A, mvalpha::metric_delsq_ck, type = type) |>
          sum(na.rm = TRUE) |>
          (\(x) x / length(K))()
    }

    denominator <- length(C) + length(K)

    if(lhs_numerator == 0 & rhs_numerator == 0 & denominator == 0){result <- 0}
    else{
      if((rlang::is_empty(C)) | rlang::is_empty(K)){result <- 1}
      else{result <- (lhs_numerator + rhs_numerator) / denominator}
    }

    return(result)

  }



#' @export

metric_delsq_ck <-
  function(c, k, type = "nominal"){as.numeric(!(c == k))}

