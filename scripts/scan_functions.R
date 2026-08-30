#########
# Functions
#########

# ── Running mean (edge-aware, O(n), NA-safe) ──────────────────────────────────
running_mean <- function(x, h) {
  n  <- length(x); if (n == 0L) return(numeric(0))
  ok <- !is.na(x)
  xc <- replace(x, !ok, 0)
  cs <- c(0, cumsum(xc)); cn <- c(0, cumsum(as.numeric(ok)))
  lo <- pmax(1L, seq_len(n) - h); hi <- pmin(n, seq_len(n) + h)
  tot <- cs[hi+1L] - cs[lo]; cnt <- cn[hi+1L] - cn[lo]
  ifelse(cnt > 0L, tot/cnt, NA_real_)
}

# ── Gap filler (mean-anchored interpolation) ─────────────────────────────────
# For each contiguous run of NAs ("gap"), compute the mean of up to h valid
# positions on each flank, then linearly interpolate across the gap between
# those two mean anchors.
#
# Why not just use the values right at the gap boundary?  The haplotype
# estimator resolves founders by cutting a distance tree (hclust + cutree).
# Positions right at the gap edge just barely passed the cutoff — their
# frequency estimates are only marginally better than the unresolved ones
# inside the gap.  Averaging h flanking positions gives robust anchor values
# driven by the well-resolved interior, not the noisy boundary.
#
# Leading/trailing NAs (no flank on one side) get a flat fill from the
# available flank's mean.  If no valid data exists at all, NAs remain.

fill_gaps <- function(x, h) {
  n  <- length(x)
  if (n == 0L || !anyNA(x)) return(x)

  ok  <- !is.na(x)
  if (!any(ok)) return(x)                   # all NA — nothing to anchor on

  # Identify contiguous NA runs (gaps)
  rle_na  <- rle(!ok)
  ends    <- cumsum(rle_na$lengths)
  starts  <- ends - rle_na$lengths + 1L

  for (g in which(rle_na$values)) {
    ga <- starts[g]                          # first NA in this gap
    gb <- ends[g]                            # last  NA in this gap

    # Left flank: mean of up to h valid positions before the gap
    left_idx <- which(ok & seq_along(x) < ga)
    if (length(left_idx) > 0L) {
      anchor_L <- mean(x[tail(left_idx, h)])
    } else {
      anchor_L <- NULL
    }

    # Right flank: mean of up to h valid positions after the gap
    right_idx <- which(ok & seq_along(x) > gb)
    if (length(right_idx) > 0L) {
      anchor_R <- mean(x[head(right_idx, h)])
    } else {
      anchor_R <- NULL
    }

    # Fill the gap
    gap_len <- gb - ga + 1L
    if (!is.null(anchor_L) && !is.null(anchor_R)) {
      # Both flanks: linear interpolation between the two trend-based anchors
      x[ga:gb] <- seq(anchor_L, anchor_R, length.out = gap_len)
    } else if (!is.null(anchor_L)) {
      # Leading edge only: flat fill from left trend
      x[ga:gb] <- anchor_L
    } else if (!is.null(anchor_R)) {
      # Trailing edge only: flat fill from right trend
      x[ga:gb] <- anchor_R
    }
    # else: no flanks at all, leave as NA
  }
  x
}


average_variance <- function(cov_matrix, tolerance = 1e-10) {
  n <- nrow(cov_matrix)  
  # Calculate eigenvalues
  eigenvalues <- eigen(cov_matrix, only.values = TRUE)$values  
  # Filter out eigenvalues that are effectively zero or negative
  positive_eigenvalues <- eigenvalues[eigenvalues > tolerance]  
  # Calculate the product of positive eigenvalues
  log_det <- sum(log(positive_eigenvalues))  
  # Use the number of positive eigenvalues for the root
  n_positive <- length(positive_eigenvalues)  
  # Calculate log of average variance
  log_avg_var <- log_det / n_positive  
  # Convert back to original scale
  avg_var <- exp(log_avg_var)  
  return(list(avg_var = avg_var, n_positive = n_positive, n_total = n))
}

wald.test3 = function(p1,p2,covar1,covar2,nrepl=1,N1=NA,N2=NA){
    
    # Wald test for multinomial frequencies
    # if nrepl = 1: (one replicate, analogous to chi square):
    #  p1 and p2 are vectors of relative frequencies to be compared
    # covar1 and covar2 are the reconstruction error 
    # covariance matrices from limSolve
    # the sampling covariance matrices are generated within limSolve
    # if nrepl > 1 (multiple replicates, analogous to CMH):
    #   p1 and p2 are matrices, each row is frequency vector for one replicate
    # covar1 and covar2 are tensors (3-dimensional arrays, third dimension 
    #  denotes replicate) for the linSolve covariance matrices
    # N1 (initial) and N2 (after treatment) 
    # are sample sizes, they are vectors when there is more than one replicate
    # N1[i], N2[i] are then for replicate i
    if (nrepl>=1){
      N1.eff=rep(NA,nrepl)
      N2.eff=rep(NA,nrepl)
      lp1 = length(p1[1,])
      cv1=array(NA,c(lp1,lp1,nrepl))
      cv2=array(NA,c(lp1,lp1,nrepl))
      for (i in 1:nrepl){
          
        covmat1  = mn.covmat((N1[i]*p1[i,]+N2[i]*p2[i,])/(N1[i]+N2[i]),2*N1[i])
        covmat2  = mn.covmat((N1[i]*p1[i,]+N2[i]*p2[i,])/(N1[i]+N2[i]),2*N2[i])
    
        N1.eff[i] = sum(diag(covmat1))*4*N1[i]^2/(sum(diag(covmat1))*2*N1[i]+2*N1[i]*sum(diag(covar1[,,i])) )
        N2.eff[i] = sum(diag(covmat2))*4*N2[i]^2/(sum(diag(covmat2))*2*N2[i]+2*N2[i]*sum(diag(covar2[,,i])) )
        cv1[,,i]= (covmat1 + covar1[ , ,i])  * (N1.eff[i])^2
        cv2[,,i]= (covmat2 + covar2[ , ,i])  * (N2.eff[i])^2
        
      }
        
      p1 = N1.eff %*% p1 / sum(N1.eff)
      p2 = N2.eff %*% p2 / sum(N2.eff)
     
      covar1= rowSums(cv1, dims = 2) / sum(N1.eff)^2
      covar2= rowSums(cv2, dims = 2) / sum(N2.eff)^2
     # browser()
    }
    else {
      covmat1  = mn.covmat((N1*p1+N2*p2)/(N1+N2),2*N1)
      covmat2  = mn.covmat((N1*p1+N2*p2)/(N1+N2),2*N2)
      covar1 = covar1 + covmat1
      covar2 = covar2 + covmat2
    }
  
  df = length(p1)-1
  covar=covar1+covar2
  eg <- eigen(covar)
  # remove last eigenvector which corresponds to eigenvalue zero
  ev <- eg$vectors[,1:df]
  eval <- eg$values[1:df]
  trafo<-diag(1/sqrt(eval)) %*% t(ev) 
  # set extremely small values to zero
  #new.covar[new.covar < 10^-9]=0
  p1= as.vector(p1); p2=as.vector(p2)
  tstat <- sum((trafo %*% (p1 - p2))^2)
  pval<- exp(pchisq(tstat,df,lower.tail=FALSE,log.p=TRUE))
  list(wald.test=tstat, p.value=pval, avg.var=average_variance(covar)$avg_var)
}

# ── Heritability from replicate haplotype frequency shifts ───────────────────
# The estimator the pipeline reports (XQTL2 #40).  Two things distinguish it
# from the older Heritability() below, which it supersedes.
#
# 1. Replicates are averaged BEFORE squaring.  They are independent measurements
#    of one shift, and mean(d^2) = mean(d)^2 + var(d), so squaring per replicate
#    and averaging after leaves the replicate scatter inside the estimate.  At
#    chrX:10,230,000 in AGE_SY20_M_no89 the old form reported 1.790 where the
#    squared mean shift is 0.505; a founder with no signal at all (t = -0.33)
#    but the largest variance at the window contributed +0.378 of it.  var(d)
#    scales as 1/N, so the spurious term roughly doubles on the male X.
#
#    Averaging first also lets the correction be MEASURED from the replicates,
#    var(d)/n, rather than modelled from mn.covmat plus the lsei covariance.
#    That matters: at that window the observed replicate scatter is 1.5-3.4x the
#    modelled variance for the three founders carrying the h2, and 0.14x for a
#    founder at C = 0.003.  A chi-square built from the measured variance gives
#    -log10p = 5.41 there against a Wald of about 5, so the test and the effect
#    size finally come from one variance rather than two.
#
# 2. The leading constant is 100*k, not 200.  For founder haplotypes at
#    frequencies p_f with additive effects a_f, truncation selection at
#    intensity i gives Dp_f = p_f * i * (a_f - abar), hence
#
#      V_A = k * sum_f p_f (a_f - abar)^2 = (k / i^2) * sum_f Dp_f^2 / p_f
#
#    with k allele copies per fly.  A hardcoded 200 asserts k = 2, which is
#    right for an autosome or a female X and wrong for a hemizygous male X
#    (k = 1) or a mixed-sex X (k = 1.5).  ploidy is k/2, i.e. xfactor, which
#    smooth_haps.R records.  Note the derivation also gives 1/p, not
#    1/(p(1-p)) -- the pq of the biallelic textbook form is what
#    p_f (a_f - abar)^2 collapses to for two alleles and does not generalise.
#
# Replicates differ in Proportion, so what they measure in common is the
# response per unit selection intensity, d/i, not d.
#
# p1, p2: nrepl x nF control and selected founder frequencies, rows in
# rep_labels order.  Returns H2 and H2_vc as percentages.
heritability_rep <- function(p1, p2, rep_labels, ProportionSelect, af_cutoff,
                             ploidy = 1.0) {
  props <- ProportionSelect$Proportion[match(rep_labels, ProportionSelect$REP)]
  ok    <- !is.na(props)
  if (sum(ok) < 2L) return(list(H2 = NA_real_, H2_vc = NA_real_))
  p1 <- p1[ok, , drop = FALSE]; p2 <- p2[ok, , drop = FALSE]
  props <- props[ok]; n <- nrow(p1)

  i_r <- dnorm(qnorm(1 - props)) / props        # selection intensity per replicate
  e   <- (p2 - p1) / i_r                        # response per unit intensity
  Cb  <- colMeans(p1)
  den <- pmax(Cb, af_cutoff)
  u   <- colMeans(e) / sqrt(den)
  s   <- (apply(e, 2, var) / n) / den

  U <- matrix(u, nrow = 1L); S <- matrix(s, nrow = 1L)
  list(H2    = 200 * ploidy * sum(u^2 - s),
       H2_vc = 200 * ploidy * vc_sum(U, S, fit_tau2(U, S)))
}

# ── Per-window variance component ────────────────────────────────────────────
# h2 is a variance, so estimate it as one.  On u = d/sqrt(p) the statistic is
# sum(u^2) and each u_hat carries known noise s, so model u ~ N(0, tau2) and fit
# tau2 by ML from the founders AT THAT WINDOW:  u_hat_f ~ N(0, tau2 + s_f).
#
# Fitting per window, not genome-wide, is the point.  With one global prior the
# ~95% of windows that are null drag tau2 to nothing and real QTL shrink about
# fivefold (0.109 against a true 0.507 in simulation, vs 0.512 for the local
# fit).  Non-negativity is then a boundary solution -- tau2 = 0 is where the
# likelihood peaks at a null window -- rather than a clamp, and founders are
# weighted by 1/(tau2+s), so a badly determined founder counts for less.
#
# Fisher-scoring fixed point, vectorised over rows.  At convergence
# sum(w^2 u^2) = sum(w) with w = 1/(tau2+s), the ML score equation.
fit_tau2 <- function(U, S, iter = 50L) {
  t2 <- pmax(rowMeans(U^2 - S), 0)
  for (k in seq_len(iter)) {
    w  <- 1 / (t2 + S)
    t2 <- pmax(rowSums(w^2 * (U^2 - S)) / rowSums(w^2), 0)
  }
  t2
}

# Posterior E[sum u^2 | u_hat] under the fitted tau2.
vc_sum <- function(U, S, t2) {
  k <- t2 / (t2 + S)
  rowSums(k^2 * U^2 + k * S)
}

mn.covmat= function(p,n,min.p=0){
  # generate multinomial covariance matrix
  # p is vector of multinomial relative frequencies
  # n is sample size
  # compute covariance matrix for relative frequencies, for absolute frequencies multiply by n^2
  # if min.p >0, then values of p smaller than min.p are set to min.p and the resulting vector is rescaled.
  p <- as.vector(p)
  p[p<min.p] = min.p; p=p/sum(p)
  mat = - tcrossprod(p)
  diag(mat) = p*(1-p)
  mat = mat/n
  mat
}
    
pseudoN.test = function(p1,p2,covar1,covar2,nrepl,N1,N2){
	pseudoN_C = rep(NA,nrepl)
	pseudoN_Z = rep(NA,nrepl)
	for(i in 1:nrepl){
		pseudoN_C[i] = (2 * N1[i] * sum(p1[i,] * (1-p1[i]))) / (2 * N1[i] * sum(diag(covar1[,,i])) + sum(p1[i,] * (1-p1[i])))
		pseudoN_Z[i] = (2 * N2[i] * sum(p2[i,] * (1-p2[i]))) / (2 * N2[i] * sum(diag(covar2[,,i])) + sum(p2[i,] * (1-p2[i])))
		}
	Count1 = round(p1*pseudoN_C,0)
	Count2 = round(p2*pseudoN_Z,0)
	lowCountFounder = apply(rbind(Count1,Count2),2,sum)
	if(sum(lowCountFounder>=5)<2){
		log10p = NA
		}else{
		Count1 = Count1[,lowCountFounder >= 5]		
		Count2 = Count2[,lowCountFounder >= 5]		
		if(nrepl==1){
			out=chisq.test(rbind(Count1,Count2),correct=TRUE)
			}else{
			nF = ncol(Count1)
			tdf = data.frame(Count=c(as.numeric(t(Count1)),as.numeric(t(Count2))),
				founder=rep(1:nF,2*nrepl),
				TRT = c(rep(1,nF*nrepl),rep(2,nF*nrepl)),
				REP = c(rep(1:nrepl,each=nF),rep(1:nrepl,each=nF)))
			D.x = xtabs(Count ~ founder + TRT + REP, data = tdf)
			out = tryCatch(mantelhaen.test(D.x,correct=TRUE),
				error=function(e) NULL)
			}
		if(is.null(out)){ return(NA) }
		log10p = -log10(out$p.value)
		}
	log10p
	}
        
add_genetic = function(df){
	df$cM = rep(NA,nrow(df))
	fm=read.table(file.path(script_dir, "../helpfiles/flymap.r6.txt"),header=FALSE)
	colnames(fm)=c("chr","pos","cM")
	library(splines)
	for(chrs in c("chrX","chr2L","chr2R","chr3L","chr3R")){
		fmX = fm %>% filter(chr==chrs)
		out = ksmooth(fmX$pos,fmX$cM,kernel="normal",bandwidth=3e6)
		f_of_x = splinefun(out$x,out$y)
		temp = f_of_x(df$pos[df$chr==chrs])
		df$cM[df$chr==chrs] = temp
		}
	df
	}

# 7-point Gauss-Hermite nodes/weights, used to integrate E[Affect^2] through the
# penetrance clamp in Heritability().  For X ~ N(mu, s^2),
#   E[f(X)] = (1/sqrt(pi)) * sum_k w_k f(mu + sqrt(2) s x_k)
GH_NODES = c(0, 0.8162878828589647, -0.8162878828589647,
             1.673551628767471, -1.673551628767471,
             2.651961356835233, -2.651961356835233)
GH_WEIGHTS = c(0.8102646175568073, 0.4256072526101278, 0.4256072526101278,
               0.05451558281912703, 0.05451558281912703,
               0.0009717812450995192, 0.0009717812450995192)

Heritability = function(p1, p2, rep_labels, ProportionSelect, af_cutoff,
			covar1=NULL, covar2=NULL, N1=NULL, N2=NULL){
	# rep_labels are the replicate LABELS for the rows of p1/p2, in row order.
	# REP is a label throughout the pipeline -- arbitrary, not necessarily
	# numeric, and not necessarily a complete 1..n run: replicates get dropped.
	# This used to rebuild REP as 1:nrepl and join ProportionSelect by that
	# invented index, which silently discarded every replicate whose label was
	# not literally its position (XQTL2 #32).  Matching on the real label makes
	# dropped replicates and non-numeric labels both work.
	nF = ncol(p1)
	nrepl = length(rep_labels)
	stopifnot(nrow(p1) == nrepl, nrow(p2) == nrepl)
	props = ProportionSelect$Proportion[match(rep_labels, ProportionSelect$REP)]
	repcol = rep(rep_labels, each=nF)
	tdf = data.frame(freq=c(as.numeric(t(p1)),as.numeric(t(p2))),
		founder=rep(1:nF,2*nrepl),
		TRT = c(rep("C",nF*nrepl),rep("Z",nF*nrepl)),
		REP = c(repcol,repcol),
		Proportion = rep(props, each=nF))

	Falconer_H2 = tdf %>%
		pivot_wider(names_from = TRT, values_from = freq) %>%
		mutate(mean_diff_sq = (Z-C)^2) %>%
		mutate(mean_af_C = case_when(C <= af_cutoff ~ af_cutoff, .default = C)) %>%
		mutate(H2temp = mean_diff_sq/mean_af_C) %>%
		group_by(REP) %>%
		summarize(H2temp_sum = sum(H2temp), Proportion = first(Proportion)) %>%
		ungroup() %>%
		filter(!is.na(Proportion)) %>%
		mutate(Falcon_i = dnorm(qnorm(1-Proportion))/Proportion) %>%
		group_by(REP) %>%
		summarize(H2 = 200 * H2temp_sum / Falcon_i^2) %>%
		ungroup() %>%
		summarize(mH2 = mean(H2)) %>%
		pull(mH2)
			
	Cutler_H2 = tdf %>%
		pivot_wider(names_from = TRT, values_from = freq) %>%
		filter(!is.na(Proportion)) %>%
		mutate(Penetrance = (Z * Proportion)/C) %>%
		mutate(Penetrance = case_when(Penetrance <= Proportion/2 ~ Proportion/2,
						  Penetrance >= 2*Proportion ~ 2*Proportion,
						  .default = Penetrance)) %>% 
		mutate(Affect = qnorm(1-Proportion) - qnorm(1-Penetrance)) %>%
		mutate(marg_Va = Affect^2 * C) %>%
		group_by(REP) %>%
		mutate(H2 = 200*sum(marg_Va)) %>%
		ungroup() %>%
		summarize(mH2 = mean(H2)) %>%
		pull(mH2)

	# ── Upward bias from squaring a noisy estimate (XQTL2 #34) ───────────────
	# Both estimators square a quantity that carries sampling error, and
	#   E[x^2] = x_true^2 + Var(x)
	# so every replicate is inflated by the same offset b.  Averaging over
	# replicates cuts the variance of H2 but not b, which is why H2 has a
	# positive floor genome-wide.  b is reported, not subtracted: H2 itself is
	# left exactly as it was.
	#
	# Var of each arm's founder frequencies = lsei reconstruction error (covar)
	# + multinomial sampling of that arm's own flies.  mn.covmat is evaluated at
	# the arm's own frequencies, not the pooled null frequencies the Wald test
	# uses, because this is the variance of the estimate actually squared.
	# C and Z are separate pools, so Cov(C,Z)=0 and no cross term appears.
	Falconer_H2_bias = NA_real_
	Cutler_H2_bias   = NA_real_
	Cutler_clamp_frac = NA_real_
	if(!is.null(covar1) && !is.null(covar2) && !is.null(N1) && !is.null(N2)){
		bF = rep(NA_real_, nrepl)
		bC = rep(NA_real_, nrepl)
		n_clamped = 0L; n_total = 0L
		for(r in seq_len(nrepl)){
			P = props[r]
			if(is.na(P)) next
			Cf = p1[r,]; Zf = p2[r,]
			vC = diag(mn.covmat(Cf, 2*N1[r]) + covar1[,,r])
			vZ = diag(mn.covmat(Zf, 2*N2[r]) + covar2[,,r])

			# Falconer squares (Z-C) directly, so the bias is exact:
			#   E[(Zhat-Chat)^2] = (Z-C)^2 + Var(Z) + Var(C)
			Falcon_i = dnorm(qnorm(1-P))/P
			bF[r] = 200 * sum((vZ + vC)/pmax(Cf, af_cutoff)) / Falcon_i^2

			# Cutler squares Affect, a nonlinear function of the frequencies:
			#   Pen    = P*Z/C          clamped to [P/2, 2P]
			#   Affect = qnorm(1-P) - qnorm(1-Pen)
			#
			# The delta method is unusable here.  Var(Pen) carries a 1/C^4
			# term, and C is at the lsei floor (3e-4) for at least one founder
			# in most windows, so the linearisation diverges -- it returned a
			# bias 14x larger than H2 itself on real data.  The clamp is what
			# actually keeps the estimator finite, and a linearisation cannot
			# see a clamp.  So integrate through it instead: Gauss-Hermite
			# quadrature of E[Affect^2] over Pen ~ N(Pen_raw, var_Pen), with the
			# clamp applied inside the integrand.  Affect is then bounded by
			# the clamp, so the bias is bounded too.
			Pen_raw = (Zf * P)/Cf
			Pen_hat = pmin(pmax(Pen_raw, P/2), 2*P)
			var_Pen = (P/Cf)^2 * vZ + (P*Zf/Cf^2)^2 * vC
			sd_Pen  = sqrt(pmax(var_Pen, 0))
			qP      = qnorm(1-P)

			EA2 = numeric(nF)
			for(k in seq_along(GH_NODES)){
				Pk = Pen_raw + sqrt(2)*sd_Pen*GH_NODES[k]
				n_clamped = n_clamped + GH_WEIGHTS[k]*sum(Pk < P/2 | Pk > 2*P)
				n_total   = n_total + GH_WEIGHTS[k]*nF
				Pk = pmin(pmax(Pk, P/2), 2*P)
				EA2 = EA2 + GH_WEIGHTS[k] * (qP - qnorm(1-Pk))^2
				}
			EA2 = EA2/sqrt(pi)

			A_hat = qP - qnorm(1-Pen_hat)
			bC[r] = 200 * sum(Cf * (EA2 - A_hat^2))
			}
		Falconer_H2_bias  = mean(bF, na.rm=TRUE)
		Cutler_H2_bias    = mean(bC, na.rm=TRUE)
		# Where the penetrance clamp binds the delta method is unreliable --
		# it linearises a hard nonlinearity.  Report how often that happened.
		Cutler_clamp_frac = if(n_total>0) n_clamped/n_total else NA_real_
		}

	list(Falconer_H2=Falconer_H2, Cutler_H2=Cutler_H2,
		Falconer_H2_bias=Falconer_H2_bias, Cutler_H2_bias=Cutler_H2_bias,
		Cutler_clamp_frac=Cutler_clamp_frac)
	}

# UNUSED in the current pipeline: hap_scan.R computes its own N1/N2 from the
# smoothed rds, where smooth_haps.R has already applied the chrX dosage factor.
# The 0.75 below is the old fixed assumption of a pool with equal numbers of
# males and females; it is wrong for a single-sex pool.  Anything reviving this
# function must take the factor from the scan's --sex the way smooth_haps.R
# does (XQTL2 #38).
doscan = function(df,chr,Nfounders){
	sexlink = 1
	if(chr=="chrX"){ sexlink=0.75 }

	# I tested with xx2$data[[1]]
	df2 = df %>%
		unnest(cols = c(sample, Groups, Haps, Err, Names)) %>%
		left_join(recodeTable) %>%
		select(-sample) %>% mutate(sample=pool) %>% select(-pool) %>%
		left_join(Numflies, join_by(sample==pool)) %>%
		separate(sample,into=c("longTRT","REP","REPrep"),remove=FALSE) %>%
		left_join(TreatmentMapping)
	
	# only analyze data for which all founders are discernable..
	allFounders = as.numeric(df2 %>% mutate(mm = max(unlist(Groups))) %>% summarize(max(mm)))	

	ll = list(Wald_log10p = NA, Pseu_log10p = NA, Falc_H2 = NA, Cutl_H2 = NA, avg.var = NA)
	if(allFounders!=Nfounders){ return(ll) }

	##  now cases where all founders are OK
	##  now collapse any pure replicates.  This is tidy ugly.  But I feel there 
	##  is value in keeping dataframe columns as lists...
	df3 = df2 %>%
		select(-Groups) %>%
		group_by(TRT,REP) %>%	
		summarise(Err_mean = list(reduce(map(Err, ~as.matrix(.x)), `+`)/length(Err)),
			Haps_mean = list(reduce(map(Haps, ~as.vector(.x)), `+`)/length(Haps)),
			Names = list(first(Names)),
			Num_mean = sexlink*mean(Num)) %>%
		rename(Haps=Haps_mean,Num=Num_mean,Err=Err_mean)

	## these summaries of the data are pretty useful for tests
	## REP is a LABEL, not an index -- arbitrary, and replicates get dropped.
	## Align the two arms on the labels present in BOTH, so a replicate missing
	## from one arm removes that pair instead of silently shifting the arms
	## against each other.  (The old code took each arm's rows in whatever order
	## they came and only compared the two counts in a dangling expression that
	## did nothing.)
	c3 = df3 %>% filter(TRT=="C")
	z3 = df3 %>% filter(TRT=="Z")
	rep_labels = intersect(c3$REP, z3$REP)
	if(length(rep_labels)==0){ return(ll) }
	c3 = c3 %>% filter(REP %in% rep_labels) %>% arrange(match(REP, rep_labels))
	z3 = z3 %>% filter(REP %in% rep_labels) %>% arrange(match(REP, rep_labels))
	nrepl = length(rep_labels)
	p1 = c3 %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t()
	row.names(p1) <- NULL
	p2 = z3 %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t()
	row.names(p2) <- NULL
	covar1 = do.call(abind, c(c3 %>% pull(Err), along = 3))
	covar2 = do.call(abind, c(z3 %>% pull(Err), along = 3))
	N1 = c3 %>% pull(Num)
	N2 = z3 %>% pull(Num)

	wt=wald.test3(p1,p2,covar1,covar2,nrepl,N1,N2)
	Wald_log10p = -log10(wt$p.value)
	Pseu_log10p = pseudoN.test(p1,p2,covar1,covar2,nrepl,N1,N2)

	af_cutoff = 0.01     # 1% --- heritability estimators can be off for really low allele frequencies
	temp = Heritability(p1, p2, rep_labels, ProportionSelect, af_cutoff,
			covar1, covar2, N1, N2)
	Falc_H2 = temp$Falconer_H2
	Cutl_H2 = temp$Cutler_H2

	ll = list(Wald_log10p = Wald_log10p, Pseu_log10p = Pseu_log10p,
			Falc_H2 = Falc_H2, Cutl_H2 = Cutl_H2, avg.var = wt$avg.var)
	ll
	}

