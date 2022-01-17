library(sf)
library(tidyverse)
library(tmap)
library(tmaptools)
library(sp)
library(INLA)

# data ####

dib <- read.table("dib2021kd.txt", header = TRUE, sep = "\t", dec = ".", na.strings = ".")
dib_sf <- st_as_sf(dib, coords=c("X.UTM20","Y.UTM20"), crs=32720)

nc.limits <- readRDS("nc.limits")
nc.limits_sf <- st_as_sf(as.data.frame(nc.limits), coords=c(c("X", "Y")), crs=22174)
nc.limits_sf <- st_transform(nc.limits_sf, crs = 32720)

dib_Kd <- dib_sf |> 
  drop_na(Kda) 

loc_Kd  <- dib_Kd |>
  st_coordinates()


dib_hl <- dib_sf |> 
  drop_na(vida_media) 

loc_hl <- dib_hl |> 
  st_coordinates()

tm_shape(dib_Kd)+
  tm_dots(shape = 25, size=1)+
tm_shape(dib_hl)+
  tm_dots(shape = 16, col = "red", size=1)

##$$\begin{align*}
# y_1(\mathbf{s}) &= \alpha_1 + z_1(\mathbf{s}) + e_1(\mathbf{s}) \\
# y_2(\mathbf{s}) &= \alpha_2 + \lambda_1 z_1(\mathbf{s}) + z_2(\mathbf{s}) + e_2(\mathbf{s}) \\
# \end{align*}$$
  
# model ####
## mesh ####

loc.obs <- st_coordinates(dib_sf)

#boundary.loc <- SpatialPoints(as.matrix(loc.obs))

boundary <- list(
  inla.nonconvex.hull(loc.obs, 60000),
  inla.nonconvex.hull(loc.obs, 90000))

## Build boundary information:
## (fmesher supports SpatialPolygons, but this app is not (yet) intelligent enough for that.)
## Build the mesh:
mesh <- inla.mesh.2d(boundary=boundary,
                     max.edge=c(20000, 80000),
                     min.angle=c(30, 20),
                     max.n=c(48000, 16000), ## Safeguard against large meshes.
                     max.n.strict=c(128000, 128000), ## Don't build a huge mesh!
                     cutoff=200, ## Filter away adjacent points.
                     offset=c(81000, 111000))

plot(mesh)
points(loc.obs, pch=16, col="red")

## spde #### with prior specification 

spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(10000, 0.01), # P(range < rangemin) = 0.01
  prior.sigma = c(10, 0.01)) # P(sigma > sigmamax) = 0.01

## hyper #### 
hyper <- list(beta = list(prior = 'normal', param = c(0, 1)))

## formula ####

form <- y ~ -1 +
  f(s1, model = spde) + 
  intercept_hl + 
  intercept_Kd +
  f(s2, model = spde) + 
  f(s12, copy = "s1", fixed = FALSE, hyper = hyper) 

## projector matrices ####


A2 <- inla.spde.make.A(mesh, loc_Kd) 
A1 <- inla.spde.make.A(mesh, loc_hl) 

#A <- inla.spde.make.A(mesh, loc.obs)

## stacks ####

stack1 <- inla.stack(
  data = list(y = cbind(log(dib_hl$vida_media), NA)),
  A = list(A1, 1), 
  effects = list(
    list(s1 = 1:spde$n.spde),
    list(intercept_Kd = rep(1, nrow(dib_hl)))
    ),
  tag = 'kd') 

stack2 <- inla.stack(
  data = list(y = cbind(NA, log(dib_Kd$Kda))),
  A = list(A2, 1), 
  effects = list(
    list(s2 = 1:spde$n.spde, s12 = 1:spde$n.spde),
    list(intercept_hl  = rep(1, nrow(dib_Kd)))),
  tag = 'hl')

stack <- inla.stack(stack1, stack2) 


# hyper.eps <- list(hyper = list(prec = list(prior = 'pc.prec', 
#                                   param = c(10000, 0.01))))
### empirical data ####
e.sd <- c(0.45,0.75)
range <- c(70000,30000)
mvar <- c(0.2,0.8)
beta <- c(0.9,0.2)

theta.ini <- c(log(1 / e.sd^2),
               c(log(range),
                 log(sqrt(mvar))), beta)
# fitting ####

result <- inla(form, rep('gaussian', 2), 
               data = inla.stack.data(stack),
               # lincomb = lincomb.predict,
               #control.family = list(hyper.eps, hyper.eps), 
               control.predictor = list(compute = TRUE,A = inla.stack.A(stack))#,
               #control.mode = list(theta = theta.ini, restart = TRUE),
               #control.inla = list(int.strategy = 'eb'),
               #control.compute = list(return.marginals.random=TRUE#,
                                      )
               #)

# result <- inla(form, data = inla.stack.data(stack), 
#             family = c('gaussian', 'gaussian', 'poisson'), 
#             control.predictor = list(compute = TRUE, 
#                                      A = inla.stack.A(stack)), 
#             control.family = list(hfix, pprec, list())) 


summary(result)

result$summary.fitted.values


sigmas1 <- inla.tmarginal(function(x) 1/sqrt(exp(x)),result$internal.marginals.hyperpar[[3]])
sigmaepsilon <- inla.tmarginal(function(x) 1/sqrt(exp(x)),result$internal.marginals.hyperpar[[1]])
restab=sapply(result$marginals.fixed, function(x) inla.zmarginal(x,silent=TRUE))
restab=cbind(restab, "alpha" = inla.zmarginal(sigmaalpha,silent=TRUE))
restab=cbind(restab, "epsilon" = inla.zmarginal(sigmaepsilon,silent=TRUE))
data.frame(restab)

ddf <- data.frame(rbind(sigmas1, sigmaepsilon),
                  errterm = gl(2, nrow(sigmaepsilon), labels = c("alpha", "epsilon")))
ggplot(ddf, aes(x, y, linetype = errterm)) + geom_line() + xlab("yield") +
  ylab("density")

# Create grid for projection
index.pred <- inla.stack.index(stack, "kd")

prj <- inla.mesh.projector(mesh)

# Settings for plotting device
par(mfrow = n2mfrow(length(names(result$summary.ran)), asp = 1))

invisible(
  sapply(names(result$marginals.fixed), function(x) {
    
    plot(result$marginals.fixed[[x]], type = 'l', 
         xlab = x, ylab = 'Density')
    abline(v = 0)
  })
)

for (j in names(result$marginals.hyperpar)) {
  ii <- result$marginals.hyperpar[[j]][,2] > sqrt(.Machine$double.eps)
  plot(result$marginals.hyperpar[[j]], 
       type = 'l', 
       xlim = range(result$marginals.hyperpar[[j]][ii, 1]), 
       xlab = names(result$marginals.hyperpar)[j], 
       ylab = 'Density',
       main = j)
  }

# Posterior mean of 'm'
par(mfrow = n2mfrow(length(names(result$summary.ran)), asp = 1))
# Posterior mean of 'm'
for (i in names(result$summary.ran)) {
  book.plot.field(field = result$summary.random[[i]]$mean, projector = prj, main = i) 
}


# 
# 
# #### INLABRU
# 
# 
# cmp <- y ~
#   -1 + 
#   covariate + 
#   field(map = coordinates, model = spde) +
#   InterceptA + 
#   InterceptB + 
#   scaling
# 
# likA <- like(
#   "normal",
#   formula = y ~ covariate + field + InterceptA,
#   data = mydata,
#   components = cmp
# )
# 
# likB <- like(
#   "normal",
#   formula = y ~ field * scaling + InterceptB,
#   data = mydata,
#   components = cmp
# )
# 
# bru(cmp, likA, likB)
# 
# 
# 
# 
