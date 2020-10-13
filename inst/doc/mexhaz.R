### R code from vignette source 'mexhaz.Rnw'

###################################################
### code chunk number 1: mexhaz.Rnw:890-893
###################################################
options(prompt = "R> ", continue = "+  ", width = 76, useFancyQuotes = FALSE)
library("mexhaz")
library("rstpm2")


###################################################
### code chunk number 2: mexhaz.Rnw:965-971
###################################################
data("simdatn1", package="mexhaz")
data("simdatn2", package="mexhaz")
simdatn1$agec <- simdatn1$age-70
simdatn2$agec <- simdatn2$age-70
head(simdatn1,3)
head(simdatn2,3)


###################################################
### code chunk number 3: mexhaz.Rnw:994-995
###################################################
Form1 <- Surv(time=timesurv, event=vstat) ~ agec + IsexH


###################################################
### code chunk number 4: mexhaz.Rnw:999-1006
###################################################
ModWb <- mexhaz(formula=Form1, data=simdatn1, base="weibull")
ModPw <- mexhaz(formula=Form1, data= simdatn1, base="pw.cst",
                knots=c(1,2,4,6,8))
ModBs <- mexhaz(formula=Form1, data=simdatn1, base="exp.bs",
                degree=3, knots=c(1,5))
ModNs <- mexhaz(formula=Form1, data=simdatn1, base="exp.ns",
                knots=c(1,5))


###################################################
### code chunk number 5: mexhaz.Rnw:1010-1018
###################################################
res <- t(sapply(list(ModWb, ModPw, ModNs, ModBs),
                function(x) {round(rbind(x$coeff[c("agec")], x$coeff[c("IsexH")],
                             x$loglik, x$n.par, -2*x$loglik+2*x$n.par),3)}))

res <- as.data.frame(res)
rownames(res) <- c("Weibull", "Piecewise Cst", "Restricted Spline", "B-Spline")
colnames(res) <- c("agec", "IsexH", "-2*log-lik", "N Param", "AIC")
res


###################################################
### code chunk number 6: mexhaz.Rnw:1055-1060
###################################################
MyTime <- seq(0,10,le=1001)
MyData <- data.frame(agec=0, IsexH=c(0,1))

P.bs.10 <- predict(ModBs, time.pts=10, data.val=MyData)
round(P.bs.10$results,3)


###################################################
### code chunk number 7: mexhaz.Rnw:1069-1071
###################################################
P.bs0 <- predict(ModBs, time.pts=MyTime, data.val=MyData[1,])
P.bs1 <- predict(ModBs, time.pts=MyTime, data.val=MyData[2,])


###################################################
### code chunk number 8: mexhaz.Rnw:1096-1103
###################################################
plot(P.bs1, which="hazard", ylim=c(0,1.5), lwd=2.5, col="blue",
main="Mortality hazard")
lines(P.bs0, which="hazard", lwd=2.5, col="red")

plot(P.bs1, which="surv", ylim=c(0,1), lwd=2.5, col="blue",
main="Overall survival")
lines(P.bs0, which="surv", lwd=2.5, col="red")


###################################################
### code chunk number 9: HazandSurv0
###################################################
par(mfrow=c(1,2), cex.lab=1.4, cex.axis=1.4, cex.main=1.4, xaxs="i", yaxs="i")

plot(NULL, xlim=c(0,10), ylim=c(0,1.5), xlab="Time", ylab="Hazard", main="Mortality hazard")
eval(parse(text="grid(lwd=2); abline(v=c(0,10),h=c(0,1.5))"))

lines(P.bs1, which="hazard", lwd=3, col="blue")
lines(P.bs0, which="hazard", lwd=3, col="red")
legend("topright",inset=0.01, c("Men", "Women"),
lwd=3, col=c("blue","red"), lty=c(1,1), cex=1.4, bg="white")

plot(NULL, xlim=c(0,10), ylim=c(0,1), xlab="Time", ylab="Survival", main="Overall survival")
eval(parse(text="grid(lwd=2); abline(v=c(0,10),h=c(0,1))"))

lines(P.bs1, which="surv", lwd=3, col="blue")
lines(P.bs0, which="surv", lwd=3, col="red")



###################################################
### code chunk number 10: HazandSurv
###################################################
P.pw0 <- predict(ModPw, time.pts=MyTime, data.val=MyData[1,])
P.ns0 <- predict(ModNs, time.pts=MyTime, data.val=MyData[1,])

par(mfrow=c(1,2), cex.lab=1.4, cex.axis=1.4, cex.main=1.4, xaxs="i", yaxs="i")

plot(NULL, xlim=c(0,10), ylim=c(0,0.7), xlab="Time", ylab="Hazard", main="Mortality hazard")
eval(parse(text="grid(lwd=2); abline(v=c(0,10),h=c(0,0.7))"))
lines(P.bs0, which="hazard", lwd=3, conf.int=FALSE)
lines(P.ns0, which="hazard", conf.int=FALSE, lwd=3, lty.pe=2)
lines(P.pw0, which="hazard", conf.int=FALSE, lwd=3, lty.pe=6)
legend("topright",inset=0.01, cex=1,
       lty = c(1, 2, 6), lwd = 3, bg="white",
       c("B-spline", "Natural spline", "Piecewise constant"),seg.len=3)

plot(NULL, xlim=c(0,10), ylim=c(0,1), xlab="Time", ylab="Survival", main="Overall survival")
eval(parse(text="grid(lwd=2); abline(v=c(0,10),h=c(0,1))"))
lines(P.bs0, which="surv", lwd=3, lty.pe=1, main="Overall survival", conf.int=FALSE)
lines(P.ns0, which="surv", conf.int=FALSE, lwd=3, lty.pe=2)
lines(P.pw0, which="surv", conf.int=FALSE, lwd=3, lty.pe=6)



###################################################
### code chunk number 11: mexhaz.Rnw:1196-1199
###################################################
ModBs2 <- mexhaz(Surv(time=timesurv, event=vstat) ~ agec + IsexH,
                 data=simdatn2, base="exp.bs", degree=2,
                 knots=c(1,5))


###################################################
### code chunk number 12: mexhaz.Rnw:1202-1205
###################################################
ModBs2.Nph <- mexhaz(Surv(time=timesurv, event=vstat)~ agec + IsexH +
nph(IsexH), data=simdatn2, base="exp.bs",
degree=2, knots=c(1,5))


###################################################
### code chunk number 13: mexhaz.Rnw:1208-1228
###################################################

ll1 <- ModBs2$loglik
ll2 <- ModBs2.Nph$loglik
np1 <- ModBs2$n.par
np2 <- ModBs2.Nph$n.par
Test <- 1-pchisq(2*(ll2-ll1),df=np2-np1)
Test

P.Bs2N.W <- predict(ModBs2.Nph, time.pts=MyTime,
                    data.val=MyData[1,], include.gradient=TRUE)
P.Bs2N.M <- predict(ModBs2.Nph, time.pts=MyTime,
                    data.val=MyData[2,], include.gradient=TRUE)

LogHR <- log(P.Bs2N.M$results$hazard) - log(P.Bs2N.W$results$hazard)

Grad.LogHR <- P.Bs2N.M$grad.loghaz - P.Bs2N.W$grad.loghaz
Var.LogHR <- diag(Grad.LogHR%*%ModBs2.Nph$vcov%*%t(Grad.LogHR))
HR <- exp(LogHR)
HR.Inf <- exp(LogHR+qnorm(0.025)*sqrt(Var.LogHR))
HR.Sup <- exp(LogHR+qnorm(0.975)*sqrt(Var.LogHR))


###################################################
### code chunk number 14: HazTD
###################################################
par(mfrow=c(1,2),cex.lab=1.4, cex.axis=1.4, cex.main=1.4, xaxs="i", yaxs="i")

plot(NULL, xlim=c(0,10), ylim=c(0,1), xlab="Time", ylab="Hazard", main="Predicted hazard at age 70")
eval(parse(text="grid(lwd=2); abline(v=c(0,10),h=c(0,1))"))
lines(P.Bs2N.M, which="hazard", ylim=c(0,0.8), lwd=3, col="blue")
lines(P.Bs2N.W, which="hazard", lwd=3, col="red")
legend("topleft",inset=0.01, c("Men", "Women"),
       lwd=3, col=c("blue","red"), lty=c(1,1), cex=1.4, bg="white")

plot(NULL, xlim=c(0,10), ylim=c(0,2.1), xlab="Time", ylab="Hazard ratio", main="Hazard ratio Men/Women")
eval(parse(text="grid(lwd=2); abline(v=c(0,10),h=c(0,2.1))"))
lines(P.Bs2N.M$results$time.pts, HR, lwd=3, type="l")
lines(P.Bs2N.M$results$time.pts, HR.Inf, type="l", lwd=2, lty="dashed")
lines(P.Bs2N.M$results$time.pts, HR.Sup, type="l", lwd=2, lty="dashed")



###################################################
### code chunk number 15: mexhaz.Rnw:1294-1297
###################################################
ModBsExc <- mexhaz(Surv(timesurv, vstat) ~ agec + IsexH +
nph(IsexH), data=simdatn1, base="exp.bs",
degree=3, knots=c(1,5), expected="popmrate")


###################################################
### code chunk number 16: mexhaz.Rnw:1300-1301
###################################################
summary(ModBsExc)


###################################################
### code chunk number 17: mexhaz.Rnw:1318-1323
###################################################
data("colon", package="rstpm2")
data("popmort", package="rstpm2")

head(colon, 3)
head(popmort, 3)


###################################################
### code chunk number 18: mexhaz.Rnw:1355-1363
###################################################
colon2 <- within(rstpm2::colon, {
    status <- ifelse(surv_mm > 120.5, 1, status)
    tm <- pmin(surv_mm, 120.5)/12
    exit <- dx + tm*365.25
    sex <- as.numeric(sex)
    exitage <- pmin(floor(age + tm), 99)
    exityear <- floor(yydx + tm)
})


###################################################
### code chunk number 19: mexhaz.Rnw:1369-1374
###################################################
colon2 <- merge(colon2, popmort, by.x=c("sex", "exitage", "exityear"),
                by.y=c("sex", "age", "year"))

head(colon2[, c("sex", "age", "stage", "status",
                "exitage", "exityear", "tm", "rate")], 3)


###################################################
### code chunk number 20: mexhaz.Rnw:1398-1410
###################################################
ModBsExc.1n <- mexhaz(Surv(timesurv, vstat) ~ agec + IsexH + nph(IsexH),
                     data=simdatn1, base="exp.bs", degree=3, knots=c(1,5),
                     n.aghq=1, expected="popmrate", random="clust")

ModBsExc.5n <- mexhaz(Surv(timesurv, vstat) ~ agec + IsexH + nph(IsexH),
                     data=simdatn1, base="exp.bs", degree=3, knots=c(1,5),
                     n.aghq=5, expected="popmrate", random="clust")

ModBsExc.10n <- mexhaz(Surv(timesurv, vstat) ~ agec + IsexH + nph(IsexH),
                     data=simdatn1, base="exp.bs", degree=3, knots=c(1,5),
                     n.aghq=10, expected="popmrate", random="clust")



###################################################
### code chunk number 21: mexhaz.Rnw:1416-1417
###################################################
summary(ModBsExc.10n)


###################################################
### code chunk number 22: mexhaz.Rnw:1425-1434
###################################################
res <- t(sapply(list(ModBsExc.1n, ModBsExc.5n, ModBsExc.10n),
                function(x) {
                  round( rbind(x$coeff[c("agec")] , exp(x$coeff[c("clust [log(sd)]")]),
                        x$loglik, -2 * x$loglik +2 * x$n.par,
                        x$time.elapsed, x$code), 3)} ))
res <- as.data.frame(res)
rownames(res) <- c("1Quad-P", "5Quad-P", "10Quad-P")
colnames(res) <- c("agec","clust (sd)","log-lik", "AIC", "time (sec)", "nlm code")
res


###################################################
### code chunk number 23: mexhaz.Rnw:1453-1455
###################################################
PBsExcR.c15 <- predict(ModBsExc.10n, time.pts=MyTime,
data.val=MyData[2,], cluster="15")


###################################################
### code chunk number 24: mexhaz.Rnw:1460-1462
###################################################
PBsExcR.0 <- predict(ModBsExc.10n, time.pts=MyTime,
data.val=MyData[2,])


