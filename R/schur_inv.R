#' The Schur method to get a inverse of a non-symmetric matrix
#'
#' @param A The matrix
#' @param extra_output logical, if you want the factorization of your matrix
#'
#' @return
#' @export
#'
#' @examples
#' schur_inv(diag(5))
schur_inv <- function(A, extra_output = FALSE){

  if(det(A) == 0) stop("\nErro!\nA matriz input é singular !\nNão faça isso!\nReveja seus conceitos")

  n_row <- dim(A)[1]
  n_col <- dim(A)[2]
  th <- ifelse(dim(A)[1]%/%2, 2, 3)

  E_m <- A[1:(n_row - th), 1:(n_col - th)]
  F_m <- A[1:(n_row - th), (n_col - th + 1):n_col]
  P_m <- A[(n_row - th + 1):n_row, 1:(n_col - th)]
  H_m <- A[(n_row - th + 1):n_row, (n_col - th + 1):n_col]

  H_1 <- solve(H_m)
  E_1 <- solve(E_m)

  S <- H_m - (P_m %*% E_1 %*% F_m)
  S_1 <- solve(S)
  T_m <- E_m - (F_m %*% H_1 %*% P_m)
  T_m_1 <- solve(T_m)


  A_1 <- matrix(0, ncol = n_col, nrow = n_row)
  A_1[1:(n_row - th), 1:(n_col - th)] <- T_m_1
  A_1[1:(n_row - th), (n_col - th + 1):n_col] <- - E_1 %*% F_m %*% S_1
  A_1[(n_row - th + 1):n_row, 1:(n_col - th)] <- - H_1 %*% P_m %*% T_m_1
  A_1[(n_row - th + 1):n_row, (n_col - th + 1):n_col] <- S_1


  if(extra_output == TRUE){
    out <- list(
      sua_matriz = A,
      suas_ABCD = list(
        T_m_1 = T_m_1,
        A_12 = - E_1 %*% F_m %*% S_1,
        A21 = - H_1 %*% P_m %*% T_m_1,
        S_1 = S_1
      ),
      sua_inversa = A_1
    )
  } else {
    out <- A_1
  }

  return(out)
}
