library(data.table)
dat <- fread("canopy-leaf-simplified.csv", header=TRUE)

# Plot BRF/albedo as function of BRF
min.wl <- 700
max.wl <- 795
with(dat[aviris_wl > min.wl & aviris_wl < max.wl],{
         x <- aviris_35_bl
         y <- x/aspen
         fit <- lm(y ~ x)
         slope <- fit$coefficients[1]
         int <- fit$coefficients[2]
         print(slope/(1-int))
         y <- x/balsampoplar
         fit <- lm(y ~ x)
         slope <- fit$coefficients[1]
         int <- fit$coefficients[2]
         print(slope/(1-int))
         y <- x/blackspruce
         fit <- lm(y ~ x)
         slope <- fit$coefficients[1]
         int <- fit$coefficients[2]
         print(slope/(1-int))
})

# Plot 
