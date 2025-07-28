set_ops <-
  function (A, B) {

    # https://stackoverflow.com/a/72631719

    ind <- match(B, A, nomatch = 0)
    A_intersect_B <- A[ind]
    A_diff_B <- A[-c(ind, length(A) + 1)]  # https://stackoverflow.com/a/52772380
    B_diff_A <- B[which(!B %in% A_intersect_B)]

    return(
      list(
        A_intersect_B = A_intersect_B,
        A_diff_B = A_diff_B,
        B_diff_A = B_diff_A
      )
    )
  }


metric_delta_CK <-
  function(C, K, type = "nominal"){

    tmp <- set_ops(C, K)

    lhs_numerator <-
      vapply(C, FUN = function(c_ind){
        vapply(tmp$B_diff_A, FUN = function(k){
          metric_delsq_ck(c_ind, k, type = type)
        }, 1) |> sum(na.rm = TRUE) # sum over k
      }, 1) |> sum(na.rm = TRUE) |> # sum over c
      (\(x) x / length(C))()

    rhs_numerator <-
      vapply(tmp$A_diff_B, FUN = function(c_ind){
        vapply(K, FUN = function(k){
          metric_delsq_ck(c_ind, k, type = type)
        }, 1) |> sum(na.rm = TRUE) # sum over k
      }, 1) |> sum(na.rm = TRUE) |> # sum over c
      (\(x) x / length(K))()

    denominator <- length(C) + length(K)

    if(length(C) == 0 | purrr::is_empty(C)) lhs_numerator <- 0
    if(length(K) == 0 | purrr::is_empty(K)) rhs_numerator <- 0
    if((length(C) == 0 | purrr::is_empty(C)) & (length(K) == 0 | purrr::is_empty(K))) denominator <- 0

    if(lhs_numerator == 0 & rhs_numerator == 0 & denominator == 0){
      result <- 0
    }else{
      if((length(C) == 0 | purrr::is_empty(C)) | (length(K) == 0 | purrr::is_empty(K))){
        result <- 1
      }else{
        result <- (lhs_numerator + rhs_numerator) / denominator
      }
    }

    return(result)

  }


metric_delsq_ck <-
  function(c, k, type = "nominal"){
    as.numeric(!(c == k))
  }

# Single Label

# cDo <-
#   vapply(units, FUN = function(u){
#     vapply(values, FUN = function(c_ind){
#       values_by_unit[, u][c_ind] *
#       vapply(values, FUN = function(k){
#         values_by_unit[, u][k] * metric_delsq_ck(c_ind, k)
#       }, 1) |> sum(na.rm = TRUE) # sum over k
#     }, 1) |> sum(na.rm = TRUE) |> # sum over c
#       (\(x) x / (m[u] - 1))()
#   }, 1) |> sum(na.rm = TRUE) |> # sum over u
#     (\(x) x / n..)()
#
# cDe <-
#      vapply(values, FUN = function(c_ind){
#       n[c_ind] *
#       vapply(values, FUN = function(k){
#         n[k] * metric_delsq_ck(c_ind, k)
#       }, 1) |> sum(na.rm = TRUE) # sum over k
#     }, 1) |> sum(na.rm = TRUE) |> # sum over c
#     (\(x) x / (n.. * (n.. - 1)))()
