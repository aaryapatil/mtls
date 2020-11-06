# Required Libraries
library('lomb')
library("multitaper")
library("LSspec")

# Artificial signal with frequencies 3 and 4 sampled at regular times
t <- seq(1, 10, length.out = 1000)
sig <- sin(2*pi*3*t) + sin(2*pi*4*t)

# Add noise to signal
noisy.sig <- sig + rnorm(length(sig), mean=0, sd=0.5)

# Plot signal
plot(t, sig, type='l', xlab = 'time', ylab = 'Signal', main = 'Artificial signal')
plot(t, noisy.sig, type='l', xlab = 'time', ylab = 'Signal', main = 'Artificial signal w noise')

# Multitaper Spectral Estimate of signal - even sampling
sigMT <- spec.mtm(sig, nw=2.0, k=5, deltat = 10.0/1000, Ftest = TRUE)
sigMTnoise <- spec.mtm(noisy.sig, nw=4.0, k=5, deltat = 10.0/1000, Ftest = TRUE)

# Plot spectrum and Ftest
plot(sigMT$freq, sigMT$spec, xlab = 'Frequency', ylab = 'P',
     main = 'Multitaper Spectral Estimate (even-sampling)', type='l')
plot(sigMT$freq, sigMT$mtm$Ftest, type = 'l', main = 'MTM F test',
     xlab = 'Frequency', ylab = 'F statistic')

# Plot spectrum and Ftest zoomed
plot(sigMTnoise$freq, sigMTnoise$spec, xlab = 'Frequency', ylab = 'P', xlim = c(0, 20),
     main = 'Multitaper Spectral Estimate (even-sampling)', type='l')
plot(sigMTnoise$freq, sigMTnoise$mtm$Ftest, type = 'l', main = 'MTM F test',
     xlab = 'Frequency', ylab = 'F statistic', xlim = c(0, 20))

# Multitaper LS Periodogram of signal
sigMTLS <- LSspecMT(t, sig, w = 4.0/1000, k = 7, Ftest=TRUE)
sigMTLSnoise <- LSspecMT(t, noisy.sig, w = 4.0/1000, k = 7, Ftest=TRUE)

# Plot spectrum and Ftest
plot(sigMTLS, type = 'l', main = 'MTLS')
sigMTLSFtest <- attr(sigMTLS, 'Ftest')
plot(sigMTLS$freq, sigMTLSFtest$F_, type = 'l', main = 'MTLS F test',
     xlab = 'Frequency', ylab = 'F statistic')

# Plot spectrum and Ftest zoomed
plot(sigMTLS, type = 'l', main = 'MTLS', xlim = c(0, 20))
plot(sigMTLS$freq, sigMTLSFtest$F_, type = 'l', main = 'MTLS F test',
     xlab = 'Frequency', ylab = 'F statistic', xlim = c(0, 20))

sigMTLSnoiseFtest <- attr(sigMTLSnoise, 'Ftest')
plot(sigMTLSnoise, type = 'l', main = 'MTLS')
# Plot Ftest
plot(sigMTLSnoise$freq, sigMTLSnoiseFtest$F_, type = 'l', main = 'MTLS F test',
     xlab = 'Frequency', ylab = 'F statistic')

# Non-uniform discrete fourier transform
sigNUMT <- NUspecMT(t, sig, w = 10.0/1000, k = 19, Ftest = TRUE)

# Plot NUDFT periodogram and Ftest
plot(sigNUMT, type = 'l', main = 'NUMT Periodogram')
# Plot Ftest
sigNUMTFtest <- attr(sigNUMT, 'Ftest')
plot(sigNUMT$freq, sigNUMTFtest$F_, type = 'l', main = 'NUMT F test',
     xlab = 'Frequency', ylab = 'F statistic')

# Plot NUDFT periodogram
plot(sigNUMT, type = 'l', main = 'NUMT Periodogram', xlim = c(0, 20))
# Plot Ftest
plot(sigNUMT$freq, sigNUMTFtest$F_, type = 'l', main = 'NUMT F test',
     xlab = 'Frequency', ylab = 'F statistic', xlim = c(0, 20))

# Willamete dataset
data(willamette)
t <- seq(1, 100, length.out = length(willamette))/2
plot(t, willamette, type='l', main = 'Willamete data', ylab = 'Signal')

# Multitaper spectral estimate
willMT <- spec.mtm(willamette, nw=4.0, k=8, deltat=100.0/length(willamette), Ftest = TRUE)
plot(willMT$freq, willMT$mtm$Ftest, type = 'l', main = 'MTM F test',
     xlab = 'Frequency', ylab = 'F statistic')

# MTLS
willMTLS <- LSspec::LSspec(t, willamette)
plot(willMTLS, type = 'l', main = 'MTLS Periodogram', xlim=c(0, 4))
willMTLSFtest <- attr(willMTLS, 'Ftest')
plot(willMTLS$freq, willMTLSFtest$F_, type = 'l', main = 'MTLS F test',
     xlab = 'Frequency', ylab = 'F statistic')

# NUDFT
willNMT <- NUspecMT(t, willamette, w = 4.0/length(willamette), k = 5, Ftest = TRUE)
plot(willNMT, type = 'l', main = 'NUMT Periodogram')
willNMTFtest <- attr(willNMT, 'Ftest')
plot(willNMT$freq, willNMTFtest$F_, type = 'l', main = 'NUMT F test',
     xlab = 'Frequency', ylab = 'F statistic')