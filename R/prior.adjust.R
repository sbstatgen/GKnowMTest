#library(locfdr)
#library(glmnet)
#library(mixtools)
#library("distr")
#library("msm")
#library("splines")


ins <- function(a, pos=c()) {
	len <- length(pos)
	for(i in 1:len) {
		nw.a <- c(a[seq(pos[i] - 1)], 1, a[seq(pos[i],length(a))])
		a <- nw.a
	}
	a
}

lfdr.alt <- function(zz1,plot=FALSE) {
	pr.del1 <- rep(NA, length(zz1))
	fsub <- is.finite(zz1)
	zz2 <- zz1[fsub]
	#pv2 <- pv[fsub]
	lres <- try(locfdr(zz2, nulltype=0, plot=plot), silent=TRUE)
		
	if(inherits(lres, "try-error"))
	{
		print("Error in locfdr. Retry with different parameters using control.locfdr or try with fit='mixnorm'")
		stop(lres)
	}
	# lres$fdr is pi0 of each snp
	pr.del1[fsub] <- (1 - lres$fdr)
	  		
	s.di <- sum(pr.del1[fsub],na.rm=TRUE)
	mu.g <- sum((pr.del1[fsub] * zz2), na.rm=TRUE) / s.di
	sig.g  <- sqrt(sum((pr.del1[fsub] * (zz2 - mu.g)^2) , na.rm=TRUE) / s.di)
	if(FALSE) {
	if(mu.g < 0) {
		lres$fdr[which(zz2 < -1)] <- 1
		pr.del1[fsub] <- (1 - lres$fdr)
	  		
		s.di <- sum(pr.del1[fsub],na.rm=TRUE)
		mu.g <- sum((pr.del1[fsub] * zz2), na.rm=TRUE) / s.di
		sig.g  <- sqrt(sum((pr.del1[fsub] * (zz2 - mu.g)^2) , na.rm=TRUE) / s.di)
	}
	}
	#sig.g <- 1
	cdf.lf <- function(x, ...) { pnorm(x,  mean=mu.g, sd=sig.g, ...) }
	pdf.lf <- function(x, ...) { dnorm(x,  mean=mu.g, sd=sig.g, ...) }
	lf = list(PPA=pr.del1, f1=pdf.lf, F1=cdf.lf)  
	lf
}

em.lfdr <- function(zz1, plot = 0)
{
	pr.del1 <- rep(NA, length(zz1))
	fsub <- is.finite(zz1)
	zz2 <- zz1[fsub]
	#pv2 <- pv[fsub]
	res <- try(locfdr(zz2, nulltype=0, plot=plot), silent=TRUE)
		
	if(inherits(res, "try-error"))
	{
		print("Error in locfdr. Retry with different parameters using control.locfdr or try with fit='mixnorm'")
		stop(res)
	}
	# lres$fdr is pi0 of each snp
	pr.del1[fsub] <- (1 - res$fdr)
	  
	#res <- locfdr(z11, nulltype=0, plot=0)
	pi0 <- res$fp0["thest", "p0"] 
	if(pi0 >= 1) pi0 <-  res$fp0["mlest", "p0"]
	if(pi0 >= 1) pi0 <-  res$fp0["cmest", "p0"]
	fdr0 <- res$fdr

	ptm1 <- proc.time()

	#zz = zz2; fdr1 = 1-fdr0 ; pi1=1-pi0; bre=120; df1=5; type="ns"; plot=0
	lfdr.alt <- function(zz, fdr1, pi1,  bre = 120, df1=5, type="ns", plot = 0)
	{
		call = match.call()

		del10 <- fdr1
		mu10 <- max(sum(zz * del10)/sum(del10), 0.0001)
		tau10 <- sqrt(max(sum((zz - mu10)^2 * del10)/sum(del10), 1.00001) - 1)
		mu.F10 <- function(xx, ...) {pnorm(xx, mean=mu10, sd=tau10, ...)}
		postr.mu <- (zz * tau10^2 + mu10)/(tau10^2 + 1)
		postr.sd <- sqrt(tau10^2/(tau10^2 + 1))

		mu.lo <- 0.0001 ; mu.up <- max(zz) + 1
		mu.breaks <- seq(mu.lo, mu.up, length.out=bre)
		# Get expected bin-counts of mu. Assuming N(mu, tau^2) as prior. 
		mu.y1 <- rep(NA, length(mu.breaks) - 1)
		for(i in 1:(length(mu.breaks)-1))
		{
			#y1[i] <- sum(fdr1[(zz > breaks[i]) & (zz <= breaks[i+1])])
			mu.y1[i] <- sum(fdr1 * (pnorm(mu.breaks[i+1], mean=postr.mu, sd=postr.sd) -  pnorm(mu.breaks[i], mean=postr.mu, sd=postr.sd)))
		}
		print(range(mu.y1))
		mu.x <- (mu.breaks[-1] + mu.breaks[ - length(mu.breaks)])/2
		mu.x <- c(mu.lo, mu.x, mu.up) 
		mu.y1 <- c(0, mu.y1, 0)
		plot(mu.x, mu.y1, type="l")
		#if(plot==1) plot(x, y1)
		if(type == "ns") {
			dat1 <- data.frame(1, ns(mu.x, df = df1))
		}
		else {
			dat1 <- data.frame(1, poly(mu.x, degree = df1))
		}
			names(dat1) <- paste("ns", 0:df1, sep="")
		fmla <- paste("mu.y1 ~0 +", paste(colnames(dat1), collapse="+"), sep="")
		res1 <- suppressWarnings(glm(formula(fmla), data=dat1, family=poisson()))
		f1 <- res1$fit
	#if(plot==1)plot(mu.x,  f1, type="l")
	
		f1.tmp <- approxfun(x=c(0, mu.x), y=c(0, f1), rule=2, ties=mean)
		f1.area <- integrate(f1.tmp, mu.lo, mu.up)$value
		fden1 <- function(xx) ifelse((xx < mu.lo) | (xx > mu.up), 0, f1.tmp(xx) / f1.area)

		if(plot==1) plot(fden1, -5, 15)
		mu.f1dist <- AbscontDistribution(d=fden1)
		mu.Fden1 <- p(mu.f1dist)
		z.f1dist <- mu.f1dist + Norm(0, 1)
		z.fden1 <- d(z.f1dist)
		z.Fden1 <- p(z.f1dist)
	
		sig10 <- sqrt(1+tau10^2)
		myF1 <- function(xx) {pnorm(xx, mean=mu10, sd=sig10)}
		if(plot==1){ 
			plot(z.Fden1, -5, 15)
			curve(myF1,  from=-5, to=15, col=3, add=TRUE)
		}
		vl = list(f1=z.fden1, F1=z.Fden1, F1.norm=function(x, ...) pnorm(x, mean=mu10, sd=sig10, ...))  
		vl$call = call
		vl
	
	}

	res1 <- lfdr.alt(zz=zz2, fdr1=1-fdr0, pi1=1-pi0,  bre = 120, df1=5, type="ns", plot = 0)

	vl = list(PPA=pr.del1, f1=res1$f1, F1=res1$F1, F1.norm=res1$F1.norm)  
    	
    	vl
}
#stop()	

RW.prior <- function(z.sc, map.obj, type="mu")
{
	agg <- aggregate(z.sc,by=list(map.obj@snp.eq), FUN=function(x) { ss <- is.finite(x); mean(x[ss])})
	z.mn <- agg[,2]
	agg <- aggregate(z.sc,by=list(map.obj@snp.eq), FUN=function(x) { ss <- is.finite(x); var(x[ss])})
	s2 <- agg[,2]
	agg <- aggregate(z.sc, by=list(map.obj@snp.eq), FUN=function(x) { ss <- is.finite(x); sum(ss)})
	nn <- agg[,2]
	pi.p <- pmin(pmax(0, ((z.mn^2)/(s2 + z.mn^2 - 1))), 1)
	mu.p <- ifelse(pi.p < (1/nn), 0 , pmax(0,z.mn) / pi.p)
	list(pp=pi.p, nn=nn, mm=mu.p,omit=NULL)

}

get.design <- function(map.obj, constant=TRUE)
{
	n.maps <- (map.obj@dim)$ob_dim["used_ob"]
	C <- (map.obj@dim)$ob_dim["no_eq_class"]
	n.anno <- (map.obj@dim)$ob_dim["no_of_anno"]

	if(n.maps > 1) {
		x.mat0 <- map.obj@eq.map
		x.mat <- matrix(0,C,n.anno)
		for(i in 1:C-1)
		{
			x.mat[i,] <- c(map.obj@eq.mat[[1]][x.mat0[i,1],], map.obj@eq.mat[[2]][x.mat0[i,2],], map.obj@eq.mat[[3]][x.mat0[i,2],])
		}
			
	} else {
		x.mat <- map.obj@eq.mat[[1]]
	}
	
	x.mat
}

#postr=pr.del1; alpha=control.prior$alpha; lambda=control.prior$lambda; use.cv=control.prior$use.cv; shrink=0
reg.prior <- function(map.obj, postr, alpha=0, lambda= NULL, shrink=0, use.cv=FALSE)
{
	n.maps <- (map.obj@dim)$ob_dim["used_ob"]
	C <- (map.obj@dim)$ob_dim["no_eq_class"]
	# pp is pi1 of each eq.class
	agg <- aggregate(postr, by=list(map.obj@snp.eq), FUN=function(x) { ss <- is.finite(x); mean(x[ss],na.rm=TRUE)})
	yy <- agg[, 2]
	agg <- aggregate(postr, by=list(map.obj@snp.eq), FUN=function(x) { ss <- is.finite(x); sum(ss,na.rm=TRUE)})
	nn <- agg[, 2]
	#yy <- postr
	pp <- rep(NA, C)
	X.mat <- get.design(map.obj, constant=TRUE)
	omit <- (!is.finite(yy))
	nn[omit] <- NA
	#use.cv <- FALSE
	if(!use.cv)
	{
		res <- glmnet(x=X.mat[!omit, ], y=cbind(1 - yy[!omit], yy[!omit]), weights=nn[!omit], offset=NULL, lambda=lambda, alpha=alpha, family="binomial")
		#res <- glm(yy[!omit] ~ X.mat[!omit, ], weights=nn[!omit])
		if(shrink == 1) pp <- predict(res, s=max(res$lambda), type="response", newx=X.mat)
		else pp[!omit] <- predict(res, s=2.312016e-05, type="response", newx=X.mat[!omit,])
	} else 
	{
		#stop("CV glmnet takes too long")
		res <- cv.glmnet(x=X.mat[!omit, ], y=cbind(1 - yy[!omit], yy[!omit]), offset=NULL, lambda=lambda, alpha=alpha, family="binomial")
		#pp.opt <- predict(res, s="lambda.min", type="response", newx=X.mat)
		pp <- predict(res, s="lambda.min", type="response", newx=X.mat)
		#pp <- predict(res, s="lambda.1se", type="response", newx=X.mat)

	}
	#pp <- pp[!is.na(pp)]
	#nn <- nn[!is.na(nn)]
	omit.pos <- which(omit==TRUE)
	if(use.cv) {
		beta.res <- coef.cv.glmnet(res,s="lambda.min")
		pp.nn <- list(pp=pp, nn=nn, omit=omit.pos, beta.res=beta.res)
	} else {
		pp.nn <- list(pp=pp, nn=nn, omit=omit.pos)
	}
	pp.nn
}


power.func1  <- function(c.vec, pp.vec, nn.vec, F1)
{
	pow <- sum(nn.vec * pp.vec * F1(c.vec, lower.tail=FALSE),na.rm=TRUE)
	pow
}

opt.pwr <- function(c.vec, pp.vec, nn.vec, F1)
{
	
	level <- 0.05
	len <- length(c.vec)
	den.w <- 1
	e.beta <- beta <- numeric(len-1)
	for(i in 1: (len-1)) {
		den <- c.vec[1]
		beta[i] <- log(c.vec[i+1]/den)
		e.beta[i] <- exp(beta[i])
		den.w <- sum(den.w, e.beta[i])
	}
	w.vec <- c(1,e.beta)
	w.vec1 <- w.vec/den.w
	c.vec1 <- qnorm((w.vec1 * level)/len, lower.tail=FALSE)
	print(nn.vec)
	print(pp.vec)
	power.func1(c.vec1, pp.vec, nn.vec, F1)

}

sim.wts <- function(prior)
{
	nn <- prior$nn
	pp <- prior$pp
	omit <- prior$omit
	if(length(omit) > 0) {
		pp <- pp[-omit]
		nn <- nn[-omit]
	}
	M <- sum(nn, na.rm=TRUE)
	#pp1 <- (1 - pp)
	if(FALSE) {
	wts <- (pp / sum(pp,na.rm=TRUE)) * M
	wts1 <- wts / nn
	}
	pi.bar <- sum((pp * nn), na.rm=TRUE) / M
	wts1 <- pp / pi.bar
	
	print(sum(wts1 * nn, na.rm=TRUE) / sum(nn, na.rm=TRUE))
	wts1
	 	
}

sim1.wts <- function(prior)
{
	nn <- prior$nn
	pp <- prior$pp
	omit <- prior$omit
	if(length(omit) > 0) {
		pp <- pp[-omit]
		nn <- nn[-omit]
	}
	M <- sum(nn, na.rm=TRUE)
	pp1 <- (1 - pp)
	wts <- (pp1 / sum(pp1,na.rm=TRUE)) * M
	wts1 <- wts / nn
	#print(mean(wts1))
	wts1
	 	
}


quad.wts <- function(prior, pr.g, level) {
	nn <- prior$nn
	pp <- prior$pp
	omit <- prior$omit
	F1 <- prior$F1
	if(length(omit) > 0) {
		pp <- pp[-omit]
		nn <- nn[-omit]
	}
	xxs <- log(pp/pr.g)
	M <- sum(nn, na.rm=TRUE)
	pwr.opt <- function(gamma.vec) {
		#gamma1 <- qnorm(gamma.vec)
		#gamma2 <- 0
		gamma1 <- gamma.vec[1]
		gamma2 <- gamma.vec[2]
		gamma3 <- gamma.vec[3]
		#gamma4 <- gamma.vec[4]
		#gamma5 <- gamma.vec[5]
		#gamma6 <- gamma.vec[6]
		#s.vec <- nn * (1 - pp) * exp((gamma1 * xxs) + (gamma2 * (xxs)^2) + (gamma3 * (xxs)^3)+ (gamma4 * (xxs)^4)  + (gamma5 * (xxs)^5) + (gamma6 * (xxs)^6))
		#med.svec <- median(s.vec)
		#scaled <- s.vec / med.svec
		#den.s <- sum(scaled, na.rm=TRUE)
		#wts.opt <- scaled /den.s
		den.s <- sum(nn * (1 - pp) * exp((gamma1 * xxs) + (gamma2 * (xxs)^2)+ (gamma3 * (xxs)^3)), na.rm=TRUE)# + (gamma4 * (xxs)^4)  + (gamma5 * (xxs)^5) + (gamma6 * (xxs)^6)), na.rm=TRUE)  
		wts.opt <- nn * (1 - pp) * exp((gamma1 * xxs) + (gamma2 * (xxs)^2)  + (gamma3 * (xxs)^3)) / den.s#+ (gamma4 * (xxs)^4) + (gamma5 * (xxs)^5) + (gamma6 * (xxs)^6)) / den.s
		#den.s <- sum(nn * exp((gamma1 * xxs) + (gamma2 * (xxs)^2)+ (gamma3 * (xxs)^3)), na.rm=TRUE)# + (gamma4 * (xxs)^4)  + (gamma5 * (xxs)^5) + (gamma6 * (xxs)^6)), na.rm=TRUE)  
		#wts.opt <- nn * exp((gamma1 * xxs) + (gamma2 * (xxs)^2)  + (gamma3 * (xxs)^3)) / den.s#+ (gamma4 * (xxs)^4) + (gamma5 * (xxs)^5) + (gamma6 * (xxs)^6)) / den.s
		crit <- qnorm(((level * wts.opt) / nn) , lower.tail=FALSE)
		o.pwr <- sum(nn * pp * F1(crit,lower.tail=FALSE), na.rm=TRUE) / sum(nn * pp, na.rm=TRUE)
		#o.pwr <- sum(nn * pp * F1(crit,lower.tail=FALSE), na.rm=TRUE) / 
		#	(sum(nn * pp * F1(crit,lower.tail=FALSE), na.rm=TRUE) + sum(nn * (1-pp) * pnorm(crit,lower.tail=FALSE), na.rm=TRUE))
		#print(range(wts.opt))
		#print(range(crit))
		#print(o.pwr)
		o.pwr
	}
	pwr.opt1 <- optim(c(0,0,0), pwr.opt, control=list(fnscale=-1))#,method="L-BFGS-B")
	gamma <- pwr.opt1$par
	#pwr.opt1 <- optimize(pwr.opt,interval=c(0,1), maximum=TRUE)
	#gamma <- c(pwr.opt1$max, 0)
	#gamma <- c(1,0)
	#s.vec <- nn * (1 - pp) * exp((gamma[1] * xxs) + (gamma[2] * (xxs)^2) + (gamma[3] * (xxs)^3) + (gamma[4] * (xxs)^4) + (gamma[5] * (xxs)^5) + (gamma[6] * (xxs)^6))
	#med.svec <- median(s.vec, na.rm=TRUE)
	#scaled <- s.vec / med.svec
	#den.s <- sum(scaled, na.rm=TRUE)
	#wts.opt <- scaled /den.s
	den.s <- sum(nn * (1 - pp) * exp((gamma[1] * xxs) + (gamma[2] * (xxs)^2) + (gamma[3] * (xxs)^3)), na.rm=TRUE)#+ (gamma[4] * (xxs)^4) + (gamma[5] * (xxs)^5) + (gamma[6] * (xxs)^6)), na.rm=TRUE)  
	wts.opt <- nn * (1 - pp) * exp((gamma[1] * xxs) + (gamma[2] * (xxs)^2) + (gamma[3] * (xxs)^3)) / den.s# + (gamma[4] * (xxs)^4) + (gamma[5] * (xxs)^5) + (gamma[6] * (xxs)^6)) / den.s
	#den.s <- sum(nn * exp((gamma[1] * xxs) + (gamma[2] * (xxs)^2) + (gamma[3] * (xxs)^3)), na.rm=TRUE)#+ (gamma[4] * (xxs)^4) + (gamma[5] * (xxs)^5) + (gamma[6] * (xxs)^6)), na.rm=TRUE)  
	#wts.opt <- nn * exp((gamma[1] * xxs) + (gamma[2] * (xxs)^2) + (gamma[3] * (xxs)^3)) / den.s# + (gamma[4] * (xxs)^4) + (gamma[5] * (xxs)^5) + (gamma[6] * (xxs)^6)) / den.s
	#wts.opt1 <- (wts.opt * M) / nn
	wts.opt1 <- wts.opt * (sum((nn * (1 - pp)), na.rm=TRUE) / (nn * (1 - pp)))
	wts.opt1
}

quad1.wts <- function(prior, pr.g, level) {
	nn <- prior$nn
	pp <- prior$pp
	omit <- prior$omit
	F1 <- prior$F1
	if(length(omit) > 0) {
		pp <- pp[-omit]
		nn <- nn[-omit]
	}
	xxs <- pp
	M <- sum(nn, na.rm=TRUE)
	pwr.opt <- function(gamma.vec) {
		gamma1 <- gamma.vec[1]
		gamma2 <- gamma.vec[2]
		gamma3 <- gamma.vec[3]
		#gamma4 <- gamma.vec[4]
		#gamma5 <- gamma.vec[5]
		#gamma6 <- gamma.vec[6]
		
		den.s <- mean(exp((gamma1 * xxs) + (gamma2 * (xxs)^2)+ (gamma3 * (xxs)^3)), na.rm=TRUE)# + (gamma4 * (xxs)^4)  + (gamma5 * (xxs)^5) + (gamma6 * (xxs)^6)), na.rm=TRUE)  
		wts.opt <- exp((gamma1 * xxs) + (gamma2 * (xxs)^2)  + (gamma3 * (xxs)^3)) / den.s#+ (gamma4 * (xxs)^4) + (gamma5 * (xxs)^5) + (gamma6 * (xxs)^6)) / den.s
		#den.s <- sum(nn * exp((gamma1 * xxs) + (gamma2 * (xxs)^2)+ (gamma3 * (xxs)^3)), na.rm=TRUE)# + (gamma4 * (xxs)^4)  + (gamma5 * (xxs)^5) + (gamma6 * (xxs)^6)), na.rm=TRUE)  
		#wts.opt <- nn * exp((gamma1 * xxs) + (gamma2 * (xxs)^2)  + (gamma3 * (xxs)^3)) / den.s#+ (gamma4 * (xxs)^4) + (gamma5 * (xxs)^5) + (gamma6 * (xxs)^6)) / den.s
		crit <- qnorm(((level * wts.opt) / nn) , lower.tail=FALSE)
		o.pwr <- sum(nn * pp * F1(crit,lower.tail=FALSE), na.rm=TRUE) / sum(nn * pp, na.rm=TRUE)
		
		o.pwr
	}
	pwr.opt1 <- optim(c(0,0,0), pwr.opt, control=list(fnscale=-1))#,method="L-BFGS-B")
	gamma <- pwr.opt1$par
	den.s <- mean(exp((gamma[1] * xxs) + (gamma[2] * (xxs)^2) + (gamma[3] * (xxs)^3)), na.rm=TRUE)#+ (gamma[4] * (xxs)^4) + (gamma[5] * (xxs)^5) + (gamma[6] * (xxs)^6)), na.rm=TRUE)  
	wts.opt <- exp((gamma[1] * xxs) + (gamma[2] * (xxs)^2) + (gamma[3] * (xxs)^3)) / den.s# + (gamma[4] * (xxs)^4) + (gamma[5] * (xxs)^5) + (gamma[6] * (xxs)^6)) / den.s
	#wts.opt1 <- wts.opt * (sum((nn * (1 - pp)), na.rm=TRUE) / (nn * (1 - pp)))
	wts.opt
}

quad2.wts <- function(prior, level=0.05, s.eq=rep(1:length(prior$nn), prior$nn))
{
	nn <- prior$nn
	pp <- prior$pp
	omit <- prior$omit
	F1 <- prior$F1
	if(length(omit) > 0) {
		pp <- pp[-omit]
		nn <- nn[-omit]
	}

	M <- sum(nn, na.rm=TRUE)
	C <- length(nn)
	d <- 3
	ww.parm <- function(parm) 
			{
			ww1 <- 1 + parm[1] * pp + parm[2] * pp^2 + parm[3] * pp^3 
			ww <- (ww1 * M)/sum(ww1 * nn, na.rm=TRUE)
			ww
			}
	ww.mfunc <- function(parm) { mean(rep(ww.parm(parm), nn), na.rm=TRUE) - 1 }
	pow.func.ww <- function(ww) {sum(nn * pp * prior$F1(qnorm(ww * level/M
									, lower.tail=FALSE), lower.tail=FALSE), na.rm=TRUE)/sum(nn * pp, na.rm=TRUE)}
	pow.func <- function(parm) {pow.func.ww(ww.parm(parm))}
	res <- try(optim(par=rep(0, d), pow.func, control=list(fnscale=-1)))
	if(inherits(res, "try-error")) stop(res)
	parm.opt <- res$par
	ww.opt <- ww.parm(parm.opt)
	#print(mean(ww.opt[s.eq],na.rm=TRUE))
	#print(pow.func(rep(0, d)))
	#print(pow.func(ww.opt))
	ww.opt
}

quant.adj <- function(z.sc, prior, snp.eq) {
	M <- length(z.sc)
	pp <- prior$pp
	nn <- prior$nn
	F1 <- prior$F1
	
	pi1g <- sum((pp * nn), na.rm=TRUE) / sum(nn, na.rm=TRUE)
	pi0g <- 1 - pi1g

	grid <- 1000000
	xx0 <- setdiff(seq(from=-1, to=1, length.out=grid),0)
	x.inf0 <- ((1/xx0) - sign(xx0))
	ord <- order(x.inf0)
	xx <- xx0[ord]
	x.inf <- x.inf0[ord]
	F1.TRUE <- F1(x.inf, lower.tail=FALSE)
	phi.TRUE <- pnorm(x.inf, lower.tail=FALSE)
	
	#stop()

	opt.h <- function(x, s.pi1, Hg.s) {
		x1 <- ((1/x) - sign(x))
		hj <- s.pi1 * F1(x1, lower.tail=FALSE) + (1 - s.pi1) * pnorm(x1, lower.tail=FALSE)
		opth <- abs(hj - Hg.s)
		opth
	}
	

	Hg <- Hj.inv <- numeric(M)
	for(i in 1:M) {
		if(is.na(z.sc[i])) { Hj.inv[i] <- NA;  next }
		F1.z <- F1(z.sc[i], lower.tail=FALSE)
		phi.z <- pnorm(z.sc[i], lower.tail=FALSE)
		Hg[i] <-  pi1g * F1.z + pi0g * phi.z
		if(Hg[i]==0) {Hg.inv[i] <- -Inf; next}
		if(Hg[i]==1) {Hg.inv[i] <- Inf; next}
		snp.pi1 <- pp[snp.eq[names(z.sc)[i]]]
		if(is.na(snp.pi1)) { Hj.inv[i] <- NA; next }
		Hj1 <- snp.pi1 * F1.TRUE + (1 - snp.pi1) * phi.TRUE
		pos.les <- Hj1 < Hg[i]
		pos.end <- which(pos.les==TRUE)[1]
		pos.str <- pos.end - 1
		int.str <- xx[pos.str]
		int.end <- xx[pos.end]
		end <- 0
		if(all(!pos.les)) { int.str <- xx[grid]; int.end <- end }
		if(all(pos.les)) { int.str <- -end; int.end <- xx[1] }
		Hj.inv0 <- try(optimize(opt.h, interval=c(int.str, int.end), s.pi1=snp.pi1
					, Hg.s=Hg[i]), silent=TRUE)
		if(inherits(Hj.inv0, "try-error")) { warning(paste("SNP:", names(z.sc)[i], Hj.inv0))
					; Hj.inv[i] <- NA }
		else { Hj.inv[i] <- ((1/Hj.inv0$min) - sign(Hj.inv0$min)) }
		print(i)
	}
	
	pcurl <- pnorm(Hj.inv, lower.tail=FALSE)
	zp.curl <- data.frame(Hj.inv, pcurl)
	zp.curl
}

power.func  <- function(w1, alpha, pp1, pp2, nn1, nn2, F1)
{
	c1 <- qnorm(alpha * w1/nn1, lower.tail=FALSE)
	c2 <- qnorm(alpha *(1-w1)/nn2, lower.tail=FALSE)
	pow <- nn1 * pp1 * F1(c1, lower.tail=FALSE) + nn2 * pp2 * F1(c2, lower.tail=FALSE)
	#pow <- pp1 * F1(c1, lower.tail=FALSE) + pp2 * F1(c2, lower.tail=FALSE)
	#print(cbind(w1, alpha, pp1, pp2, pow))
	pow
}

#prior <- list(F1=function(x, ...) pnorm(q=x, mean=2, sd=3, ...), F1inv=function(y, ...) qnorm(p=y, mean=2, sd=3, ...), nn=11:1000, pp=runif(990, 0.001, 0.4))
#s.eq=rep(1:length(prior$nn), prior$nn)
#level=0.05; npth.sc=-9
#prior$pp <- rep(prior$pp[1],length(prior$pp))
pair.wts <- function(prior, level=0.05, npth.sc=0)
{
	F1 <- prior$F1
	n.vec <- prior$nn
	pi.vec <- prior$pp
	npth.wn <- which.max(n.vec)
	npth.n <- n.vec[npth.wn]
	if(npth.sc >= 0) {
		omit <- c(npth.wn, prior$omit)
	} else { 
		omit <- prior$omit
	}
	if(length(omit) > 0) {
		
		pi.vec <- pi.vec[-omit]
		n.vec <- n.vec[-omit]
	} 

	
	C <- length(pi.vec)
	#ord <- order(n.vec,  pi.vec)
	ord <- order(pi.vec)

	M <- sum(n.vec,na.rm=TRUE)
	if(npth.sc >= 0) {
	M <- M + npth.n 
	npth.alpha <- 10^(-npth.sc) * level
	alpha.cur <- (level - npth.alpha)
	} else {
	alpha.cur <- level
	}
	#eps <- 0.0000001
	eps0 <- 1e-30
	alpha.vec <- rep(NA, C)
	crit.vec <- rep(NA, C)
	for( i in 1:(C-1))
	{
		sub1 <- ord[i]
		sub2 <- ord[(i + 1):C]
		nn1 <- sum(n.vec[sub1],na.rm=TRUE)
		nn2 <- sum(n.vec[sub2],na.rm=TRUE)
		pp1 <- sum(n.vec[sub1]*pi.vec[sub1],na.rm=TRUE)/nn1
		pp2 <- sum(n.vec[sub2]*pi.vec[sub2],na.rm=TRUE)/nn2
		eps <- (eps0 * nn1)
		res <- optimize(function(w1) power.func(w1, alpha.cur, pp1, pp2, nn1, nn2, F1)
				, interval=c(eps, 1 - eps), maximum=TRUE, tol=1e-10)
		w1 <- res$maximum
		pow1 <- power.func(w1, alpha.cur, pp1, pp2, nn1, nn2, F1)
		pow2 <- power.func((nn1/(nn1 + nn2)), alpha.cur, pp1, pp2, nn1, nn2, F1)
		cat(c(i, pp1, pp2, nn1, nn2, w1, nn1/(nn1 + nn2), pow1, pow2), "\n", append=(i > 1), file="debug.txt")
		alpha.vec[ord[i]] <- w1 * alpha.cur
		crit.vec[ord[i]] <- qnorm(alpha.vec[ord[i]]/nn1, lower.tail=FALSE)
		alpha.cur <- (1 - w1) * alpha.cur
	}
	alpha.vec[ord[C]] <- alpha.cur
	crit.vec[ord[C]] <- qnorm(alpha.vec[ord[C]]/n.vec[ord[C]], lower.tail=FALSE)
	if(npth.sc > 0) {
		alpha.vec <- c(alpha.vec, npth.alpha)
		crit.npth <- qnorm(npth.alpha/npth.n, lower.tail=FALSE)
		crit.vec <- c(crit.vec, crit.npth)
		n.vec <- c(n.vec, npth.n)
	}
	wt.vec <- (M / level) * pnorm(crit.vec, lower.tail=FALSE)
	#print(sum(wt.vec * n.vec,na.rm=TRUE)/sum(n.vec,na.rm=TRUE))
	wt.vec1 <- wt.vec / (sum(wt.vec * n.vec,na.rm=TRUE)/sum(n.vec,na.rm=TRUE))
	print(sum(wt.vec1 * n.vec,na.rm=TRUE)/sum(n.vec,na.rm=TRUE))
	#print(wt.vec[length(wt.vec)])
	wt.vec1
}

		
opt.wts <- function(prior, level=0.05)
{
	nn <- prior$nn
	mm <- prior$mm
	pp <- prior$pp
	omit <- prior$omit
	if(length(omit) > 0) {
		pp <- pp[-omit]
		nn <- nn[-omit]
	}
	pp1 <- pp
	pp0 <- (1 - pp)
	M <- sum(nn, na.rm=TRUE)
	F1 <- prior$F1
	get.lam1 <- function(lam1, mm, pp, nn, level)
	{
		lam <- 10^lam1
		crit <- (((mm^2) / 2) - log(pp / lam)) / (mm)
		type1 <- sum(nn * pnorm(crit, lower.tail=FALSE),na.rm=TRUE)
		lamd <- (type1 - level)
		return(lamd)
	}

	lam1.res <- try(uniroot(get.lam1, interval=c(-1000, 1000),tol=1e-10, level=level, mm = mm, pp = pp, nn = nn),silent=TRUE)
	print(lam1.res)
	lam <- 10^lam1.res$root
	if(inherits(lam1.res, "try-error")) {
		stop("err")
		grid <- 1000000
		#lam1 <- seq(length.out=grid,0,1)
		#lam.inf <- (1/lam1) - 1
		lam1 <- seq(length.out=grid,lam1.res$root-0.15,lam1.res$root+0.15)
		lam.inf <- 10^lam1
		lam.alpha <- numeric(grid)
		lamd <- numeric(grid)
		for(i in 1:grid) {
			crit <- (((mm^2) / 2) - log(pp / lam.inf[i])) / (mm)
			lam.alpha[i] <- sum(nn * pnorm(crit, lower.tail=FALSE))
			lamd[i] <- (if(lam.alpha[i] < level) 1 else 0) * (lam.alpha[i])	
		}
		lam <- lam.inf[which.max(lamd)]
	}

	
	c.opt <- (((mm^2) / 2) - log(pp / lam)) / (mm)
        wt.opt <- (nn * pnorm(c.opt, lower.tail=FALSE))/sum(nn * pnorm(c.opt, lower.tail=FALSE),na.rm=TRUE)
	wt.opt1 <- wt.opt * (M/nn)
	print(sum(wt.opt1 * nn,na.rm=TRUE)/sum(nn,na.rm=TRUE))	
	wt.opt1
}

opt1.wts <- function(prior, level=0.05, s.eq=rep(1:length(prior$nn), prior$nn))
{
	nn <- prior$nn
	mm <- prior$mm
	pp <- prior$pp
	omit <- prior$omit
	if(length(omit) > 0) {
		pp <- pp[-omit]
		nn <- nn[-omit]
	}
	M <- sum(nn, na.rm=TRUE)
	C <- length(nn)
	#ww.c <- function(cc) {(M/level) * pnorm(prior$F1inv(pmin(1, cc * pp), lower.tail=FALSE), lower.tail=FALSE)}
	ww.c <- function(cc) {(M/level) * pnorm((mm)/2 + log(cc/(mm * pp)),lower.tail=FALSE)}
	ww.mfunc <- function(cc) { mean(rep(ww.c(cc), nn), na.rm=TRUE) - 1 }
	pow.func.ww <- function(ww) {sum(nn * pp * prior$F1(qnorm(ww * level/M, lower.tail=FALSE)
							, lower.tail=FALSE), na.rm=TRUE)/sum(nn * pp, na.rm=TRUE)}
	pow.func.cc <- function(cc) {pow.func.ww(ww.c(cc))}
	res <- try(uniroot(ww.mfunc, interval=c(0, 10000), tol=1e-10, check.conv=TRUE))
	if(inherits(res, "try-error")) stop(res)
	cc.opt <- res$root
	ww.opt <- ww.c(cc.opt)
	avg <- mean(ww.opt[s.eq], na.rm=TRUE)
	print(avg)
	if(abs(avg - 1)> 0.001) stop("Mean weight is not 1 !")
	#print(pow.func(rep(0, d)))
	#print(pow.func(ww.opt))
	ww.opt
}


#stop()
#prior.meth="Reg";adj.meth="pair.wt";control.prior=list(alpha=1,  lambda=NULL, alt.dist="locfdr", use.cv=FALSE)
#control.adj=list(level=0.05); verbose=TRUE; control.locfdr=list(plot=FALSE); control.pair=list(npth.sc=-1)
prior.adjust <- function(data.df, map.obj, prior.meth=c("MoM","Reg"), adj.meth=c("quantile","opt.wt","simple.wt","quad.wt","pair.wt")
		, control.prior=list(alpha=1,  lambda=NULL, alt.dist="norm", use.cv=FALSE), control.adj=list(level=0.05)
		, verbose=TRUE, control.locfdr=list(plot=FALSE), control.pair=list(npth.sc=-1)) {
	# getting the data #	
	if("zscore" %in% names(data.df)) {
		zz <- data.df$zscore 
		pv <- 2 * pnorm(abs(zz), lower.tail=FALSE)
		zz1 <- qnorm(pv, lower.tail=FALSE)
	} else if("pvalue" %in% names(data.df)) {
		pv <- data.df$pvalue
		zz1 <- qnorm(pv, lower.tail=FALSE)
	} else { 
		zz <- data.df[, "Beta"]/data.df[, "SD"]
		pv <- 2 * pnorm(abs(zz), lower.tail=FALSE)
		zz1 <- qnorm(pv, lower.tail=FALSE)
	}
	ptm0 <- proc.time()
	names(zz1) <- rownames(data.df)

	fsub <- is.finite(zz1)
	zz2 <- zz1[fsub]

	# model fitting #
	if(verbose) print("Fitting the global model ...")
	if(control.prior$alt.dist=="norm") {

	   	alt.dist <- lfdr.alt(zz1,plot=control.locfdr$plot)
		pr.del1 <- alt.dist$PPA
		prior <- list(F1=alt.dist$F1,f1=alt.dist$f1)

	} else if(control.prior$alt.dist=="npar") {
		
		alt.dist <- em.lfdr(zz1)
		pr.del1 <- alt.dist$PPA
		prior <- list(F1=alt.dist$F1,f1=alt.dist$f1 )
		#curve(alt.dist$F1(x),-10,10)
		#curve(alt.dist$F1.norm(x),-10,10, add=TRUE, col=2)
		

	} else if(control.prior$alt.dist=="mixnorm") {
		
		fsub <- is.finite(zz1)
		zz2 <- zz1[fsub]
		out <- try(normalmixEM(zz2, mean.constr=c(0, NA), sd.constr=c(1, NA)), silent=TRUE)
		if(inherits(out, "try-error"))
		{
			print("Error in normalmixEM. Retry with different parameters using control. locfdr or try with alt.distF='mixnorm'")
			stop(out)
		}
		mu.g <- out$mu[2]
		sig.g <- min(1, out$sigma[2])
		pi.1 <- out$lambda[2]
		pi.0 <- (1 - pi.1)
		cdf.mxn <- function(x, ...) { pnorm(x,  mean=mu.g, sd=sig.g, ...) }
		pdf.mxn <- function(x, ...) { dnorm(x,  mean=mu.g, sd=sig.g, ...) }
		pr.del1 <- 1/(1 +  exp(log(pi.0) - log(pi.1) + dnorm(zz, log=TRUE) - dnorm(zz, mean=mu.g, sd=sig.g, log=TRUE)))
		prior <- list(F1=cdf.mxn,f1=pdf.mxn)
	
	} else {
		stop("alt.dist method should be either 'locfdr' or 'mixnorm'")
	}
	
	# get posterior from input #
	if("PPA" %in% colnames(data.df)) { pr.del1 <- data.df$PPA }

	mu.g <- sum(zz2 * pr.del1[fsub], na.rm=TRUE)/sum(pr.del1, na.rm=TRUE)
	if(FALSE) {
	if(mu.g < 0) {
		lres$fdr[which(zz2 < -1)] <- 1
		pr.del1[fsub] <- (1 - lres$fdr)
	  		
		s.di <- sum(pr.del1[fsub],na.rm=TRUE)
		mu.g <- sum((pr.del1[fsub] * zz2), na.rm=TRUE) / s.di
	}
	}
	pr.g <- mean(pr.del1[fsub], na.rm=TRUE)
	
#stop()
	
	# obtaining prior #
	if(verbose) print("Obtaining Prior ...")
	if(prior.meth=="MoM") {
		res <- RW.prior(z.sc=zz1, map.obj)
		prior$mm <- res$mm
		#pr.g <- mean(pr.del1, na.rm=TRUE)
	} else if(prior.meth=="Reg") {
		res <- reg.prior(map.obj,  postr=pr.del1, alpha=control.prior$alpha, lambda=control.prior$lambda, use.cv=control.prior$use.cv)
		#z.mn <- aggregate(zz1,by=list(map.obj@snp.eq), FUN=function(x) { ss <- is.finite(x); mean(x[ss], na.rm=TRUE)})
		#mu.p <- z.mn$x / res$pp
		#pr.g <- mean(pr.del1, na.rm=TRUE) # for locfdr
		#mu.g <- mean(is.finite(zz1),na.rm=TRUE) / pr.g
		#pr.del1.1 <- pr.del1[fsub]
		#mu.g <- sum(zz2 * pr.del1.1)/sum(pr.del1.1) # for locfdr
		#beta.res <- res$beta.res
		prior$mm <- mu.g
	} else {
		stop("prior method should be either 'method of moments(MoM)' or 'Regression(Reg)'")
	}
	prior <- c(prior, list(pp=res$pp), list(nn=res$nn), list(omit=res$omit))
#stop()
	# getting adjusted weights #
	if(adj.meth=="quantile") {
		names(map.obj@snp.eq) <- names(zz1)
		Padj <- quant.adj(z.sc=zz1, prior=prior, snp.eq=map.obj@snp.eq)
		M <- nrow(map.obj@snp.df)
		summary.mat <- matrix(NA, M, 4)
		rownames(summary.mat) <- rownames(map.obj@snp.df)
		colnames(summary.mat) <- c("Z", "P", "Z.adj", "P.adj")
		summary.mat[, "Z"] <- zz1
		summary.mat[, "P"] <- pv
		summary.mat[, "Z.adj"] <- Padj[,1]
		summary.mat[, "P.adj"] <- Padj[,2]
	} else {
		if(adj.meth=="opt.wt") { wts <- opt1.wts(prior, level=control.adj$level) }
		if(adj.meth=="simple.wt") { wts <- sim1.wts(prior) }
		if(adj.meth=="quad.wt") { wts <- quad.wts(prior, pr.g, control.adj$level) }
		#if(adj.meth=="quad.wt") { wts <- quad2.wts(prior, control.adj$level) }
		if(adj.meth=="pair.wt") { wts <- pair.wts(prior, level=control.adj$level, npth.sc=control.pair$npth.sc) }
		if(length(res$omit) > 0) {
			wts <- ins(wts, pos=res$omit)
		}
		ptm1 <- proc.time()	
		if(verbose) { print("Completed weighted analysis"); print(ptm1 - ptm0) }
		M <- nrow(map.obj@snp.df)
		PPA <- rep(NA,M)
		prior.f1 <- prior$f1
		prior1 <- prior$pp[map.obj@snp.eq][fsub]
		prior0 <- 1 - prior1
		PPA[fsub] <- 1 / (1 + exp(-log(prior.f1(zz2)) - log(prior1) + dnorm(zz2,log=TRUE) + log(prior0)))
		summary.mat <- matrix(NA, M, 5)
		rownames(summary.mat) <- rownames(map.obj@snp.df)
		colnames(summary.mat) <- c("Z", "P", "wP", "adjP","PPA")
		summary.mat[, "Z"] <- zz1
		summary.mat[, "P"] <- pv
		summary.mat[, "wP"] <- pmin(summary.mat[, "P"] / wts[map.obj@snp.eq], 1)
		summary.mat[, "adjP"] <- pmin(M * summary.mat[, "P"] / wts[map.obj@snp.eq], 1)
		summary.mat[, "PPA"] <- PPA

	}
	#res.df <- list(res.df0=cbind(map.obj@snp.df,  summary.mat), beta.res=beta.res)
	#res.df <- cbind(map.obj@snp.df,  summary.mat)
	res.df <- list(res.df0=cbind(map.obj@snp.df,  summary.mat), pp = prior$pp[map.obj@snp.eq], alt.d = prior.f1, wt.vec = wts[map.obj@snp.eq])
	res.df
}

