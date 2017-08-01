set.seed(123)
beta  = c(rnorm(10,0,2.5), rep(0, 32))
pmat1sided <- matrix(1-pnorm(rnorm(40,beta)), nrow=2, ncol=20, byrow=TRUE) #one-sided p-values
pmat2sided = 2*pmin(pmat1sided, 1-pmat1sided) #two-sided p-values

#compute the one-sided Fisher global null p-value and the conditional p-values given that Fisher's
# global null p-value is at most 0.001
out.Fisher= conditional.pvalues(pmat1sided, globaltest = "Fisher", pthreshold=0.001)
print(out.Fisher)

#plot the original p-values, and the conditional p-values after global null selection by Fisher,
# for the selected row.
par(mfrow=c(1,1))
plot(seq(1:20), -log10(pmat1sided[1,]), col="red", ylab = "-log10(PV)",
     main="conditional and original one-sided p-values on -log10 scale" )
points(seq(1:20), -log10(out.Fisher$p2C[1,]), pch = 2)
legend(16,8, c("original",  "cond Fisher"), pch = c(1,2), col=c("red", "black"))
abline(-log10(0.05/20),0,lty=2, col="gray")


#compute the two-sided Fisher global null p-value and the conditional p-values given that Fisher's
# global null p-value is at most 0.001
out.Fisher= conditional.pvalues(pmat2sided, globaltest = "Fisher", pthreshold=0.001)
print(out.Fisher)

#plot the original p-values, and the conditional p-values after global null selection by Fisher,
# for the selected row.
par(mfrow=c(1,1))
plot(seq(1:20), -log10(pmat2sided[1,]), col="red", ylab = "-log10(PV)",
     main="conditional and original two-sided p-values on -log10 scale" )
points(seq(1:20), -log10(out.Fisher$p2C[1,]), pch = 2)
legend(16,8, c("original",  "cond Fisher"), pch = c(1,2), col=c("red", "black"))
abline(-log10(0.05/20),0,lty=2, col="gray")



#compute the Pearson global null p-value and the conditional p-values given that Pearson's
# global null p-value is at most 0.001
out.Pearson = conditional.pvalues(pmat1sided, globaltest = "Pearson", pthreshold=0.001)
print(out.Pearson)

#plot the original p-values, and the conditional p-values after global null selection by Pearson,
# for the selected row.
par(mfrow=c(1,1))
plot(seq(1:20), -log10(pmat2sided[1,]), col="red", ylab = "-log10(PV)",
     main="conditional and original two-sided p-values on -log10 scale" )
points(seq(1:20), -log10(out.Pearson$p2C[1,]), pch = 2)
legend(16,8, c("original",  "cond Pearson"), pch = c(1,2), col=c("red", "black"))
abline(-log10(0.05/20),0,lty=2, col="gray")

