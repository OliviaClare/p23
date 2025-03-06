#' p value combination approach using closed testing procedure (only up to 3 dose levels)
#'
#' This functions implements p-value combination approach using closed testing procedure for phase 2/3 inferential design with multiple dose levels at Stage 1.
#' Multiple testing adjustment method is Simes method for each family wise hypothesis H_J at stage 1. Then closed testing procedure is applied for weighted z combining stage 1 and stage 2.
#' 
#' 
#' @param z1 A matrix of z values at stage 1. z1 has dimensions of (ntrials, ndoses). ndoses must be 2 or 3.
#' @param p1 A matrix of p values at stage 1. p1 = 1-pnorm(z1). either z1 or p1 must be provided.
#' @param z2 stage 2 z test statistic. a vector of length = ntrials. length(z2) must be equal to nrow(z1)
#' @param p2 stage 2 p value. p2 = 1-pnorm(z2). Either z2 or p2 must be provided.
#' @param w weight for stage 1 z statistic. If w is a single number, all trials use the same weight in p value combination.
#' @param bd.z rejection boundary in z scale for the combination test
#' @param method "simes", "Dunnett". Currently only "simes" method is implemented.
#' 
#' @return Returned values include:
#' \describe{
#' \item{selected.dose}{Selected dose}
#' \item{unadj.ps}{unadjusted p value for selected dose s}
#' \item{adj.ps}{adjusted p value for selected dose s in Stage 1 after closed testing procedure using Simes method}
#' \item{comb.p}{combination p value combining stage 1 and stage 2}
#' \item{unadj.zs}{unadjusted z value for selected dose s}
#' \item{adj.zs}{adjusted z value for selected dose s in Stage 1 after closed testing procedure using Simes method}
#' \item{comb.z}{Combination z value combining stage 1 and stage 2}
#' \item{RejPct}{Percentage of rejection, applicable for simulations when z1 has a large number of rows (trials)}
#' \item{Selection}{Probability of selection of each dose, applicable for simulations when z1 has a large number of rows (trials)}
#' }
#' 
#' @examples
#' 
#' #Example (1): Stage 1 has 3 z statistics: z1 = c(0.784, 1.118, 1.941)
#' Select 3rd dose to proceed to Stage 2. At Stage 2, the incremental z test statistic is z2 = 2.1 and the rejection boundary is 1.96.
#' 
#' 
#' ph23.comb.p(z1=matrix(c(0.784, 1.118, 1.941), nrow=1),  z2 = 2.1, bd.z=1.96, w=0.2)
#' 
#' ph23.comb.p(z1=matrix(c(0.784, 1.118), nrow=1),  z2 = 2.1, bd.z=1.96, w=0.2)
#' 
#' ph23.comb.p(z1=matrix(c(0.784), nrow=1),  z2 = 2.1, bd.z=1.96, w=0.2)
#' 
#' #Example (2): check type I error by simulations
#' nSim = 100000
#' z1 = matrix(NA, nrow=nSim, ncol=3)
#' for (j in 1:3){z1[ ,j] = rnorm(nSim)}
#' z2 = rnorm(nSim)
#' 
#' o=ph23.comb.p(z1=z1,  z2 = z2, bd.z=1.96, w=0.2)
#' 
#' #Example (3). w is random between 0 and 1
#' o=ph23.comb.p(z1=z1,  z2 = z2, bd.z=1.96, w=runif(nSim))
#' 
#' @export 
#' 
ph23.comb.p = function(z1, z2, p1=NULL, p2=NULL, bd.z=1.96, w=0.2, method="simes"){
  n.doses = ncol(z1) #only works for up to 3 dose levels in this program.
  N = nrow(z1) #number of trials
  
  if (!is.null(p1)) {z1 = qnorm(1-p1)}
  if (!is.null(p2)) {z2 = qnorm(1-p2)}
  
  #familywise P(reject H_J | H_J) assuming independent incremental
  
  ##################
  #Stage 1
  ##################
  zs = rep(NA, N) #raw z for Stage 1 selected dose s 
  s = rep(NA, N) #Stage 1 selected dose s
  
  for (i in 1:N){
    tmp = sort(z1[i,], index.return = TRUE)
    zs[i] = tmp$x[n.doses]
    s[i] = tmp$ix[n.doses]
  }
  #cbind(z, zs, ps, s)
  ps = 1-pnorm(zs) #raw p value for selected dose
  
  #Simes method for adjustment
  
  p = 1 - pnorm(z1)
  simes.p3 = rep(NA, N) #adjusted p by all doses for stage 1 on the selected dose s
  if (n.doses > 1) {simes.p2 = matrix(NA, nrow=N, ncol=n.doses-1)} #adjusted p by all 2 doses that includes s
  p1s = rep(NA, N) #final adjusted p for stage 1 on the selected dose s with Closed testing procedure
  
  if (n.doses > 1) {
  for (i in 1:N){
    #adjustment for all p values
    simes.p3[i] = simes(p[i,])
    
    #adjustment for any 2 p values that includes s
    unselected = (1:n.doses)[-s[i]]
      for (j in 1:length(unselected)){
        simes.p2[i,j] = simes(c(p[i, s[i]], p[i,unselected[j]]))
      }
    
    #final adjusted p by CTP
    p1s[i] = max(simes.p3[i], simes.p2[i,], p[i, s[i]])
  }
  } else {
    p1s = ps
  }
  
  z1s = qnorm(1-p1s) #adjusted z for stage 1 on the selected dose s
  #cbind(s, p, simes.p3, simes.p2, ps, p1s, z1s)
  
  ##################
  #Stage 2, independent incremental z for the selected dose
  ##################

  #Final test statistic
  comb.z = w*z1s + sqrt(1-w^2)*z2
  comb.p = 1-pnorm(comb.z)
  
  o = list()
  o$selected.dose = s
  o$unadj.ps = ps
  o$adj.ps = p1s
  o$comb.p = comb.p
  
  o$unadj.zs = zs
  o$adj.zs = z1s
  o$comb.z = comb.z
  
  selection = rep(NA, n.doses)
  if (N > 1){ 
    RejPct = sum(comb.z > bd.z) / N
    o$RejPct = RejPct
    for (j in 1:n.doses){
      selection[j] = sum(s == j) / N
    }
    o$selection = selection
    
  }
  
  return(o)
}

