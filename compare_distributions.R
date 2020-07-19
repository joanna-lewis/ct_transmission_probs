rm(list=ls())

load("../ct_transmission_prob_results/transmission_probability_reparam_natsal2_hetnonew_1000.RData")
fit_hetnonew <- fit_all
op_hnn <- extract(fit_hetnonew)

load("../ct_transmission_prob_results/transmission_probability_reparam_natsal2_nonocon_1000.RData")
fit_nonocon <- fit_all
op_nnc <- extract(fit_nonocon)

####################
# plot
####################

h1_hnn <- hist(op_hnn$rho_mf, breaks = c(seq(0,1,0.05), max(op_hnn$rho_mf)), plot=FALSE)
h2_hnn <- hist(op_hnn$rho_fm, breaks = c(seq(0,1,0.05), max(op_hnn$rho_fm)), plot=FALSE)

h1_nnc <- hist(op_nnc$rho_mf, breaks = c(seq(0,1,0.05), max(op_nnc$rho_mf)), plot=FALSE)
h2_nnc <- hist(op_nnc$rho_fm, breaks = c(seq(0,1,0.05), max(op_nnc$rho_fm)), plot=FALSE)

#quartz()
dev.off()
par(mfrow=c(1,1), mar=c(4,4,2,2))
plot(rep(h1_hnn$breaks, each=2), c(0, rep(h1_hnn$density, each=2), 0), 
     type='l', col='#F98400', lwd = 3, 
     xlim=c(0,1), ylim=c(0,5), xaxt = "n",
     xlab="Transmission probability per partnership (%)", ylab='Density'
)
lines(rep(h2_hnn$breaks, each=2), c(0, rep(h2_hnn$density, each=2), 0), col='#00A08A', lwd=3)
axis(1, at = seq(0,1,0.2), labels=seq(0,100,20))
legend('topright', 
       lwd=3, lty=c(1,1,1,2),
       col=c('#F98400', '#00A08A', 'grey', 'grey'), 
       legend=c('Male-to-female', 'Female-to-male', 'New partners', 'Partners without\na condom'), 
       bty='n',
       cex=0.8)

lines(rep(h1_nnc$breaks, each=2), c(0, rep(h1_nnc$density, each=2), 0), col='#F98400', lwd=3, lty=2)
lines(rep(h2_nnc$breaks, each=2), c(0, rep(h2_nnc$density, each=2), 0), col='#00A08A', lwd=3, lty=2)



quantile(op_0$rho_mf, p=c(0.5, 0.025, 0.975))
quantile(op_0$rho_fm, p=c(0.5, 0.025, 0.975))
