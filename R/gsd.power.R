#' Calculate power for group sequential design by simulation
#'
#' This functions calculates the cumulative power and overall power for a group sequential design by simulations.
#'
#' @param z Standard z test statistics. A matrix of dimensions (n_trials, n_analyses).
#' @param p p values (1-sided),i.e., p = 1-pnorm(z)
#' @param bd.z z value boundary
#' @param bd.p p value boundary
#' 
#' @return Cumulative power by analysis
#' 
#' @examples
#' #Example (1): Stage 1: 4 arms; 3 dose levels; each arm 50 patients.
#' #Stage 2: additional 200 patients per arm will be enrolled at stage 2
#' #medians for the 4 arms: 9, 11, 13 and control = 8 months
#' #Enrollment: 12 months uniform in stage 1; 12 months uniform in stage 2
#' #Holding period: 4 months between stage 1 and 2
#' #Dose selection will be based on data cut at 16 months
#' #Stage 2 has 2 planned analyses at 300 and 380 events respectively.
#'
#' #Using O'Brien Fleming boundary, the rejection boundaries are: 
#' bd.z = actualBounds(planned.events=c(300, 380), act.events=c(300, 380), sf=gsDesign::sfLDOF, alpha=0.025)$actual.z
#' #2.268527 2.022098
#' 
#' o = simu.ph23data(nSim=1000, n1 = c(50, 50, 50, 50), n2 = c(200, 200), m = c(9, 11, 13, 8), 
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' DCO1 = 16, Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4, targetEvents2 = c(300, 380))
#' 
#' IA = ph23.comb.p(z1=o$z1,  z2 = o$z2[,1], bd.z=bd.z[1], w=o$w[,1])
#' FA = ph23.comb.p(z1=o$z1,  z2 = o$z2[,2], bd.z=bd.z[2], w=o$w[,2])
#' 
#' #power calculation
#' gsd.power(z = cbind(IA$comb.z, FA$comb.z), bd.z=bd.z)
#' 
#' @export 
#' 
gsd.power = function(z, bd.z, p=NULL, bd.p=NULL){
  
  if (!is.null(z)){K = ncol(z)} else {K = ncol(p)}
  if (!is.null(z)){ # filter out NA rows in z and corresponding bd.z YC ============
    narows = which(is.na(rowSums(z)))
    if(length(narows)!=0){
      warning(paste0("A total of ", length(narows), " trials have no interim analysis."))
      z = z[-narows,]
      bd.z = bd.z[-narows,]
    }
    N = nrow(z)
  } else {N = nrow(p)}
  # if (!is.null(z)){N = nrow(z)} else {N = nrow(p)}
  cum.pow = rep(NA, K)
  
  if (!is.null(z)) {
    if(is.matrix(bd.z)){
      for (j in 1:K){
        tmp = 0
        for (ii in 1:j){
          tmp = tmp + as.numeric(z[,ii] > bd.z[,ii]) # altered to reflect actual bd.z YC=========
        }
        cum.pow[j] = sum(tmp > 0) / N
      }
    }else{
      for (j in 1:K){
        tmp = 0
        for (ii in 1:j){
          tmp = tmp + as.numeric(z[,ii] > bd.z[ii]) # one boundary for all (original)
        }
        cum.pow[j] = sum(tmp > 0) / N
      }
    }
    
  } else {
    for (j in 1:K){
      tmp = 0
      for (ii in 1:j){
        tmp = tmp + as.numeric(p[,ii] < bd.p[ii])
      }
      cum.pow[j] = sum(tmp > 0) / N
    }
  }
  
  return(cum.pow)
}

