#' Using a theorem to get the inverse of a non-symmetric matrix
#'
#' @param W The input matrix
#' @param extra_output logical, if you want the factorization of your matrix
#'
#' @return
#' @export pracma::rref
#'
#' @examples
#' theorem_1_9(matrix(nrow = 5, ncol = 5, runif(5^2)))
theorem_1_9 <- function(W, extra_output = FALSE){

  if(det(W) == 0) stop("\nErro!\nA matriz input é singular !\nNão faça isso!\nReveja seus conceitos")
  Inde <- diag(nrow(W))
  if(nrow(W) != ncol(W)) stop("\nA matriz input não é quadrada!\nO teorema não lida com esse tipo de matriz.")
  if(all.equal(W, Inde)){
    if(extra_output == TRUE) warning(paste0("A matriz input é a matriz identidade n = ",nrow(W),
                                            ". As matrizes A, B, C e D não serão printadas"))
    out <- list(W)
  } else {
    #Quantas linhas são Linearmente independentes ----
    reduc_depen <- pracma::rref(W - Inde)
    n_depen <- dim(reduc_depen[rowSums(reduc_depen),])[1]

    # SVD ----
    S_svd <- svd(W - Inde,
                 # Ja que a matriz e quadrada:
                 nu = n_depen, nv = n_depen)
    S_sig <- diag(S_svd$d)
    S_sig_f <- S_sig[rowSums(S_sig) > 1e-10, colSums(S_sig) > 1e-10]
    S_v <- S_svd$v
    S_u <- S_svd$u

    #S_svd_output <- zapsmall(S_u %*% S_sig_f %*% t(S_v))

    out_inverse <- (Inde - Inde %*% S_u %*% solve(solve(S_sig_f)+ t(S_v) %*% Inde %*% S_u) %*% t(S_v) %*% Inde)

    out_inverse <- zapsmall(out_inverse)

    if(extra_output == TRUE){
      out <- list(sua_matriz = W,
                  suas_ABCD = list(
                    A = diag(nrow(W)),
                    B = S_sig_f,
                    C = S_u,
                    D = t(S_v)
                  ),
                  sua_inversa = out_inverse
      )

    } else {
      out <- out_inverse
    }
  }
  return(out)

}
