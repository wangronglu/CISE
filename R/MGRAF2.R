#' Second variant of M-GRAF model
#'
#' \code{MGRAF2} returns the estimated common structure Z and \eqn{\Lambda} that
#' are shared by all the subjects as well as the subject-specific low rank
#' matrix \eqn{Q_i} for multiple undirected graphs.
#'
#' The subject-specific deviation \eqn{D_i} is decomposed into \deqn{D_i = Q_i
#' * \Lambda * Q_i^{\top},} where each \eqn{Q_i} is a VxK orthonormal matrix and
#' \eqn{\Lambda} is a KxK diagonal matrix.
#'
#' @param A Binary array with size VxVxn storing the VxV symmetric adjacency
#'   matrices of n graphs.
#' @param K An integer that specifies the latent dimension of the graphs
#' @param tol A numeric scalar that specifies the convergence threshold of CISE
#'   algorithm. CISE iteration continues until the absolute percent change in
#'   joint log-likelihood is smaller than this value. Default is tol = 0.01.
#' @param maxit An integer that specifies the maximum number of iterations.
#'   Default is maxit = 5.
#'
#' @return A list is returned containing the ingredients below from M-GRAF2
#'   model corresponding to the largest log-likelihood over iterations.
#'   \item{Z}{A numeric vector containing the lower triangular
#'   entries in the estimated matrix Z.} \item{Lambda}{Kx1 vector storing the
#'   diagonal entries in \eqn{\Lambda}.} \item{Q}{VxKxn array containing the estimated
#'   VxK orthonormal matrix \eqn{Q_i}, i=1,...,n.} \item{D_LT}{Lxn matrix where
#'   each column stores the lower triangular entries in \eqn{D_i = Q_i * \Lambda
#'   * Q_i^{\top}}; L=V(V-1)/2.} \item{LL_max}{Maximum log-likelihood across
#'   iterations.} \item{LL}{Joint log-likelihood at each iteration.}
#'
#' @examples
#' data(A)
#' res = MGRAF2(A=A, K=5, tol=0.01, maxit=5)
#'
#' @import gdata
#' @import Matrix
#' @import glmnet
#' @import rARPACK
#' @importFrom stats coef sd
#' @export

####################################################################################################################################### CISE algorithm to estimate common structure and low-dimensional individual
####################################################################################################################################### structure of multiple undirected binary networks #####
##------- M-GRAF Model, Variant 2 -----------------------------------##
## A_i ~ Bernoulli(\Pi_i) logit(\Pi_i) = Z + D_i = Z + Q_i %*% \Lambda %*% t(Q_i)
## Lambda is the same for all subjects The algorithm iterate between the following
## steps until convergence 1. Given Z, \Lambda_i, solve for Q_i by doing
## eigen-decomposition on (A_i-P_0), where P_0 = 1/(1+exp(-Z)) 2. Given Q_i, solve Z
## and \Lambda_i by logistic regression

MGRAF2 = function(A, K, tol, maxit) {

    n = dim(A)[3]
    V = dim(A)[1]
    L = V * (V - 1)/2

    if (missing(tol)) {
        tol = 0.01
    }

    if (tol<=0){
      stop("threshold tol should be positive")
    }

    if (missing(maxit)) {
        maxit = 5
    }
    ###### Initialization ----------------------------------------------------------------

    ## initialize P_0 by A_bar
    ## -----------------------------------------------------------
    P0 = apply(A, c(1, 2), sum)/n
    ## initialize Z by log odds of P0
    vec_P0 = lowerTriangle(P0)
    vec_P0[which(vec_P0 == 1)] = 1 - (1e-16)  # note 1 - (1-(1e-17)) == 0 in R
    vec_P0[which(vec_P0 == 0)] = 1e-16
    Z = log(vec_P0/(1 - vec_P0))

    ## initialize Lambda by averaging eigenvalues of (A_i-P0) (sorted decreasingly in
    ## magnitude)
    Lambda = apply(A, 3, function(x) {
        ED = eigs_sym(A = x - P0, k = K, which = "LM", opts = list(retvec = FALSE))  # eigenvalues only
        return(ED$values)
    })  # Kxn

    if (K == 1) {
        Lambda = mean(Lambda)
    } else {
        Lambda = apply(Lambda, 1, mean)
    }

    ## initialize Q_i by eigenvectors of (A_i-P0)
    Q = array(0, c(V, K, n))
    # number of positive lambda's in Lambda
    pl = sum(Lambda > 0)

    for (i in 1:n) {
        if (pl == 0) {
            ED = eigs_sym(A = A[, , i] - P0, k = K, which = "SA")
            Q[, , i] = ED$vectors
        } else if (pl == K) {
            ED = eigs_sym(A = A[, , i] - P0, k = K, which = "LA")
            Q[, , i] = ED$vectors
        } else {
            ED1 = eigs_sym(A = A[, , i] - P0, k = pl, which = "LA")
            ED2 = eigs_sym(A = A[, , i] - P0, k = K - pl, which = "SA")
            Q[, , i] = cbind(ED1$vectors, ED2$vectors)
        }
    }

    ###---------- compute initial log-likelihood ---------------------------------------------------
    A_LT = apply(A, 3, lowerTriangle)

    ## an array of lower-triangle of principal matrices specific to subject
    M_array = apply(Q, c(2, 3), function(x) {
        lowerTriangle(tcrossprod(x))
    })  # LxKxn
    D_LT = matrix(0, nrow = L, ncol = n)

    LL_A = 0

    for (i in 1:n) {
        if (K == 1) {
            D_LT[, i] = M_array[, , i] * Lambda
        } else {
            D_LT[, i] = M_array[, , i] %*% Lambda
        }

        vec_Pi = 1/(1 + exp(-Z - D_LT[, i]))

        vec_Pi[which(vec_Pi == 1)] = 1 - (1e-16)  # note 1 - (1-(1e-17)) == 0 in R
        vec_Pi[which(vec_Pi == 0)] = 1e-16
        LL_A = LL_A + sum(A_LT[, i] * log(vec_Pi) + (1 - A_LT[, i]) * log(1 - vec_Pi))
    }

    ###################################################################################### TUNE PENALTY PARAMETER LAMBDA IN GLMNET ptm = proc.time() CONSTRUCT Y
    ###################################################################################### ----------------------------------------------------------------------
    y = factor(c(A_LT))

    ### CONSTRUCT PENALTY FACTORS FOR Z AND LAMBDA ---------------------------------------
    ### prior precision of Z
    phi_z = 0.01
    # prior precision of lambda
    s_l = 2.5  # prior scale
    phi_lambda = 1/(s_l^2)
    # penalty factor
    pen_fac = c(rep(phi_z, L), rep(phi_lambda, K))
    # normalize to ensure sum(pen_fac) = L+K, #variables
    const_pf = sum(pen_fac)/(L + K)
    pen_fac = pen_fac/const_pf
    # glmnet penalty factor
    lambda_glm = c(10^(0:-8), 0) * const_pf

    ### CONSTRUCT DESIGN-MATRIX ----------------------------------------------------------
    ### construct intercept part of design matrix
    ### -----------------------------------------
    design_int = Diagonal(L)
    for (i in 2:n) {
        design_int = rbind(design_int, Diagonal(L))
    }

    ## construct predictors M part of design matrix
    ## -------------------------------------- check s.d of vec(M^i_k)
    sd_M = apply(M_array, 2, sd)  # Kx1

    M_list = lapply(1:n, function(i) {
        # scale M_array[,k,i] to have sd 0.5
        if (K == 1) {
            temp_M = M_array[, , i]/2/sd_M  # Lx1
        } else {
            temp_M = sweep(M_array[, , i], 2, 2 * sd_M, FUN = "/")  # LxK
        }
    })

    if (K == 1) {
        design_mat = cbind(design_int, do.call(c, M_list))
    } else {
        design_mat = cbind(design_int, do.call(rbind, M_list))
    }

    rm(M_list)

    ## run cv.glmnet to determine optimal penalty lambda
    rglmModel = cv.glmnet(x = design_mat, y = y, family = "binomial", alpha = 0, lambda = lambda_glm,
        standardize = FALSE, intercept = FALSE, penalty.factor = pen_fac, maxit = 200,
        nfolds = 5, parallel = FALSE)  # type.measure=deviance
    ind_lambda_opt = which(lambda_glm == rglmModel$lambda.min)
    glm_coef = coef(rglmModel, s = "lambda.min")[-1]

    ##----- update Z and P0 -----##
    Z = glm_coef[1:L]  # Lx1
    P0 = matrix(0, nrow = V, ncol = V)
    lowerTriangle(P0) = 1/(1 + exp(-Z))
    P0 = P0 + t(P0)

    ##----- update Lambda -----##
    Lambda = glm_coef[(L + 1):(L + K)]  # Kx1
    # unscale Lambda
    Lambda = Lambda/sd_M/2
    # sort lambda
    if (K > 1) {
        Lambda = sort(Lambda, decreasing = TRUE)
    }

    #################################################### 2-step Iterative Algorithm #####################
    LL_seq = numeric(maxit + 1)
    LL_seq[1] = LL_A

    # elapse_time = numeric(maxit)

    for (st in 1:maxit) {
        ptm = proc.time()
        ########### Update Q ------------------------------------------------------------------ num of
        ########### positive lambda's in Lambda
        pl = sum(Lambda > 0)

        Q = array(0, c(V, K, n))

        for (i in 1:n) {
            if (pl == 0) {
                ED = eigs_sym(A = A[, , i] - P0, k = K, which = "SA")
                Q[, , i] = ED$vectors
            } else if (pl == K) {
                ED = eigs_sym(A = A[, , i] - P0, k = K, which = "LA")
                Q[, , i] = ED$vectors
            } else {
                ED1 = eigs_sym(A = A[, , i] - P0, k = pl, which = "LA")
                ED2 = eigs_sym(A = A[, , i] - P0, k = K - pl, which = "SA")
                Q[, , i] = cbind(ED1$vectors, ED2$vectors)
            }
        }

        ######### COMPUTE JOINT LOGLIKELIHOOD -----------------------------------------------------
        ######### an array of lower-triangle of principal matrices specific to subject
        M_array = apply(Q, c(2, 3), function(x) {
            lowerTriangle(tcrossprod(x))
        })  # LxKxn
        D_LT = matrix(0, nrow = L, ncol = n)

        LL_A = 0

        for (i in 1:n) {
            if (K == 1) {
                D_LT[, i] = M_array[, , i] * Lambda
            } else {
                D_LT[, i] = M_array[, , i] %*% Lambda
            }

            vec_Pi = 1/(1 + exp(-Z - D_LT[, i]))

            vec_Pi[which(vec_Pi == 1)] = 1 - (1e-16)  # note 1 - (1-(1e-17)) == 0 in R
            vec_Pi[which(vec_Pi == 0)] = 1e-16
            LL_A = LL_A + sum(A_LT[, i] * log(vec_Pi) + (1 - A_LT[, i]) * log(1 - vec_Pi))
        }

        LL_seq[st + 1] = LL_A

        print(st)

        if (LL_seq[st + 1] > max(LL_seq[1:st])) {
            D_LT_best = D_LT
            Q_best = Q
            Lambda_best = Lambda
            Z_best = Z
            LL_max = LL_seq[st + 1]
        }

        if (abs(LL_seq[st + 1] - LL_seq[st])/abs(LL_seq[st]) < tol) {
            break
        }

        ############ CONSTRUCT DESIGN-MATRIX FOR LOGISTIC REGRESSION ----------------------------------
        ############ intercept part of design matrix has been constructed
        ############ -------------------------------------- construct predictors M part of design
        ############ matrix -------------------------------------- scale M
        sd_M = apply(M_array, 2, sd)  # Kx1

        M_list = lapply(1:n, function(i) {
            # scale M_array[,k,i] to have sd 0.5
            if (K == 1) {
                temp_M = M_array[, , i]/2/sd_M  # Lx1
            } else {
                temp_M = sweep(M_array[, , i], 2, 2 * sd_M, FUN = "/")  # LxK
            }
        })

        if (K == 1) {
            design_mat = cbind(design_int, do.call(c, M_list))
        } else {
            design_mat = cbind(design_int, do.call(rbind, M_list))
        }

        rm(M_list)

        ########## LOGISTIC REGRESSION ------------------------------------------------------------
        ########## run a penalized logistic regression (ridge regression) Instead of setting penalty
        ########## = lambda_opt, we use a sequence of larger penalty parameters as warm starts.  This
        ########## is more robust though may take longer time.

        rglmModel = glmnet(x = design_mat, y = y, family = "binomial", alpha = 0, lambda = lambda_glm[1:ind_lambda_opt],
            standardize = FALSE, intercept = FALSE, penalty.factor = pen_fac, maxit = 200)

        ind_beta = dim(rglmModel$beta)[2]
        ##----- update Z and P0 -----##
        Z = rglmModel$beta[1:L, ind_beta]  # Lx1
        P0 = matrix(0, nrow = V, ncol = V)
        lowerTriangle(P0) = 1/(1 + exp(-Z))
        P0 = P0 + t(P0)

        ##----- update Lambda -----##
        Lambda = rglmModel$beta[(L + 1):(L + K), ind_beta]  # Kx1
        # unscale Lambda
        Lambda = Lambda/sd_M/2
        # sort lambda
        if (K > 1) {
            Lambda = sort(Lambda, decreasing = TRUE)
        }

        # elapse_time[st] = as.numeric((proc.time()-ptm))[3]
    }

    results = list(Z = Z_best, Lambda = Lambda_best, Q = Q_best, D_LT = D_LT_best, LL_max = LL_max,
        LL = LL_seq)
    return(results)
}
