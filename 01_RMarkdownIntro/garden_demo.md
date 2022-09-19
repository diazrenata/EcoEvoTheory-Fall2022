Revisiting the garden SAR
================

-   <a href="#load-data" id="toc-load-data">Load data</a>
-   <a href="#wrangling" id="toc-wrangling">wrangling</a>

Code all taken from the code for the Why Theory lecture, removing saving
plots and a couple of dependencies (upscale, socorro).

### Load data

``` r
# garden SAR ----
g <- read.csv(here::here("01_RMarkdownIntro", "data", "garden.csv"), as.is = TRUE, header = FALSE)
```

### wrangling

``` r
glist <- unlist(g)
names(glist) <- NULL 
glist <- strsplit(glist, '; ')


# site by spp matrix
allPlant <- sort(unique(unlist(glist)))


gmat <- lapply(1:length(glist), function(i) {
    as.integer(allPlant %in% glist[[i]])
})
gmat <- do.call(rbind, gmat)




m <- matrix(1:64, nrow = 8)
j <- rep(2^(0:3), c(1, 2, 2, 2))
i <- rep(2^(0:3), c(2, 2, 2, 1))

S <- sapply(1:length(j), function(k) {
    ind <- as.vector(m[1:i[k], 1:j[k]])
    sum(colSums(gmat[ind, , drop = FALSE]) > 0)
})


dat <- data.frame(area = i * j * 0.3 * 0.4, S = S)
dat$logA <- log(dat$area, 10)
dat$logS <- log(dat$S, 10)
```

``` r
# linear SAR build-up
pdf(here::here("01_RMarkdownIntro", "outputs", 'fig_SAR_build-up1.pdf'), width = 4, height = 4)

par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(1.75, 0.4, 0), tcl = -0.35, 
    cex = 1.2, lwd = 1.2, 
    bg = 'black', fg = 'white', col.axis = 'white', col.lab = 'white')
plot(dat$area, dat$S, ylim = c(0, 16), xlim = c(0, 8),
     xaxt = 'n', type = 'n', 
     xlab = 'Area', ylab = 'Number of Species', 
     cex = 1.2, lwd = 2)
axis(1, at = 10 * (0:6) * (0.3 * 0.4), labels = 10 * (0:6))

dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# linear SAR build-up
pdf(here::here("01_RMarkdownIntro", "outputs", 'fig_SAR_build-up2.pdf'), width = 4, height = 4)

par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(1.75, 0.4, 0), tcl = -0.35, 
    cex = 1.2, lwd = 1.2, 
    bg = 'black', fg = 'white', col.axis = 'white', col.lab = 'white')
plot(dat$area, dat$S, ylim = c(0, 16), xlim = c(0, 8),
     xaxt = 'n', 
     xlab = 'Area', ylab = 'Number of Species', 
     cex = 1.2, lwd = 2)
axis(1, at = 10 * (0:6) * (0.3 * 0.4), labels = 10 * (0:6))

dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# linear SAR
pdf(here::here("01_RMarkdownIntro", "outputs", 'fig_SAR_linear.pdf'), width = 4, height = 4)
par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(1.75, 0.4, 0), tcl = -0.35, 
    cex = 1.2, lwd = 1.2, 
    bg = 'black', fg = 'white', col.axis = 'white', col.lab = 'white')
plot(dat$area, dat$S, xlim = c(0, 8), ylim = c(0, 16), 
     xlab = expression('Area ('*m^2*')'), ylab = 'Number of Species', 
     cex = 1.2, lwd = 2)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# loglog SAR
pdf(here::here("01_RMarkdownIntro", "outputs", 'fig_SAR_loglog.pdf'), width = 4, height = 4)
par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(1.75, 0.4, 0), tcl = -0.35, 
    cex = 1.2, lwd = 1.2, 
    bg = 'black', fg = 'white', col.axis = 'white', col.lab = 'white')
plot(dat$area, dat$S, log = 'xy', axes = FALSE,
     xlim = c(0.08, 10), ylim = c(0.8, 30),
     xlab = expression('Area ('*m^2*')'), ylab = 'Number of Species', 
     cex = 1.2, lwd = 2)

#socorro::logAxis(1:2)
box()

#dev.off()


# loglog SAR with line
pdf(here::here("01_RMarkdownIntro", "outputs", 'fig_SAR_final.pdf'), width = 4, height = 4)
par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(1.75, 0.4, 0), tcl = -0.35, 
    cex = 1.2, lwd = 1.2, 
    bg = 'black', fg = 'white', col.axis = 'white', col.lab = 'white')
plot(dat$area, dat$S, log = 'xy', axes = FALSE,
     xlim = c(0.08, 10), ylim = c(0.8, 30),
     xlab = expression('Area ('*m^2*')'), ylab = 'Number of Species', 
     cex = 1.2, lwd = 2)

#socorro::logAxis(1:2)
box()
m <- lm(logS ~ logA, data = dat)
curve(10^m$coefficients[1] * x^m$coefficients[2], lwd = 2,
      from = 0.1, to = 10, add = TRUE)

dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# garden METE ----

library(meteR)

rc <- expand.grid(1:8, 1:8)
ii <- rep(1:nrow(rc), sapply(glist, length))
rc <- rc[ii, ]

meteDat <- data.frame(spp = unlist(glist), row = rc[, 1], col = rc[, 2])
sar <- meteSAR(S0 = 16, N0 = 250, Amin = 1, A0 = 64)

pdf(here::here("01_RMarkdownIntro", "outputs", 'fig_gardenMETE.pdf'), width = 4, height = 4)

par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(1.75, 0.4, 0), tcl = -0.35, 
    cex = 1.2, lwd = 1.2, 
    bg = 'black', fg = 'white', col.axis = 'white', col.lab = 'white')
plot(dat$area, dat$S, log = 'xy', axes = FALSE,
     xlim = c(0.08, 10), ylim = c(0.8, 30),
     xlab = expression('Area ('*m^2*')'), ylab = 'Number of Species', 
     cex = 1.2, lwd = 2)

#socorro::logAxis(1:2)
box()
m <- lm(logS ~ logA, data = dat)
curve(10^m$coefficients[1] * x^m$coefficients[2], lwd = 2,
      from = min(dat$area), to = max(dat$area), add = TRUE)

lines(sar$pred$A * 0.3 * 0.4, sar$pred$S, col = hsv(0.5, 1, 0.8), lwd = 2)

dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# faking z slopes ----

set.seed(1)
y <- runif(200, 0.1, 0.8)
x <- runif(200)

pdf(here::here("01_RMarkdownIntro", "outputs", 'fig_z_build-up.pdf'), width = 4.5, height = 4.5)
par(bg = 'black', fg = 'white', col.axis = 'white', col.lab = 'white')
plot(x, y, xlim = c(0, 8), axes = FALSE, lwd = 0.6)
box(lwd = 0.6)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# 
# 
# # upscale prediction ----
# 
# # data taken from figure in john's book
# xy <- read.csv(here::here("01_RMarkdownIntro", "outputs", 'upscale.csv', as.is = TRUE)
# xmod <- lm(y ~ x, data = xy[xy$type == 'xaxis', ])
# ymod <- lm(y ~ x, data = xy[xy$type == 'yaxis', ])
# 
# # converting pixle coords to plot coords
# upSAR <- xy[xy$type %in% c('mete', 'obs'), ]
# upSAR$x <- xmod$coefficients[1] + xmod$coefficients[2] * upSAR$x
# upSAR$y <- ymod$coefficients[1] + ymod$coefficients[2] * upSAR$y
# 
# 
# # smooth mete fun
# smoothSAR <- as.data.frame(spline(upSAR[upSAR$type == 'mete', 2:3]))
# smoothSAR <- cbind(type = 'mete', smoothSAR)
# upSAR <- rbind(smoothSAR, upSAR[upSAR$type == 'obs', ])
# 
# # transforming to linear scale
# upSAR$x <- exp(upSAR$x)
# upSAR$y <- exp(upSAR$y)
# 
# 
# 
# plot(upSAR[upSAR$type == 'mete', 2:3], col = 'red', type = 'l', log = 'xy')
# 
# points(upSAR[upSAR$type == 'obs', 2:3])
```
