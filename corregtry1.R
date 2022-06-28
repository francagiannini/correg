library(sf)
library(tidyverse)
library(tmap)
library(tmaptools)
library(sp)
library(INLA)
library(wesanderson)
library(maptiles)

# data ####

dib <- read.table("dib2021kd.txt", header = TRUE, sep = "\t", dec = ".", na.strings = ".")

dib_sf <- st_as_sf(dib, coords=c("X.UTM20","Y.UTM20"), crs=32720)

nc.limits <- readRDS("nc.limits")
nc.limits_sf <- st_as_sf(as.data.frame(nc.limits), coords=c(c("X", "Y")), crs=22174)
nc.limits_sf <- st_transform(nc.limits_sf, crs = 32720)

limits_sf <- st_read("C:\\Users\\au710823\\Dropbox\\Franca\\Doctorado\\Tesis\\codigo_resumen_tesis\\archivos\\limites_cba_cut.shp")

zones_sf <- st_read(
  "C:\\Users\\au710823\\Dropbox\\Franca\\Doctorado\\Tesis\\MapeoDigitalCordoba-master\\MapeoDigitalCordoba-master\\Datos/zonas_capI.kml")|> 
  st_transform( crs = 32720)
# 
# osmnolbl <- list(
#   src = 'osmnolabel',
#   q = 'https://{s}.tiles.wmflabs.org/osm-no-labels/{z}/{x}/{y}.png',
#   sub = c('a','b', 'c'),
#   cit = '© OpenStreetMap contributors.'
# )
# # dowload tiles and compose raster (SpatRaster)
# 
# mtqOSM <- maptiles::get_tiles(x = limits_sf, 
#                               provider = osmnolbl, 
#                               crop = TRUE,
#                               zoom = 6)
# # Plot the tiles
# tilesLayer(mtqOSM)



dib_Kd <- dib_sf |> 
  drop_na(Kda) 

loc_Kd  <- dib_Kd |>
  st_coordinates()

dib_hl <- 


dib_hl <- dib_sf |> 
  drop_na(vida_media) 

loc_hl <- dib_hl |> 
  st_coordinates()

tmap_mode("plot")

# tm_shape(limits_sf)+
#   tm_polygons(col="#F2F2F2")+
tm_shape(zones_sf)+
  tm_polygons(col="#F2F2F2")+
tm_shape(dib_Kd)+
  tm_dots(shape = 25, size=0.3, col="#FFA700")+
tm_shape(dib_hl)+
  tm_dots(shape = 16, col = "#226E53", size=0.5)

my_binwidth <- 5 

ggplot(data, aes(x)) +        # Draw hist & density with count on y-axis
  geom_histogram(binwidth = my_binwidth) +
  geom_density(aes(y = ..density.. * (nrow(data) * my_binwidth)), col = 2)

hl_his <- 
  dib_sf |> ggplot(aes(vida_media))+
  #geom_histogram(binwidth = my_binwidth) +
  geom_density(aes(y = ..density.. * (nrow(dib_hl) * my_binwidth)), fill="#226E53")+
  #theme_bw()+
  theme(panel.background = NULL)+
  labs(x="half-life", y= "Frequency")+
  geom_vline(xintercept=16.69)#mean(dib_hl$vida_media)
  
kd_his <- 
  dib_sf |> ggplot(aes(Kda))+
  #geom_histogram(binwidth = 1) +
  geom_density(aes(y = ..density.. * (nrow(dib_Kd) * 1)), fill="#FFA700")+
  #theme_bw()+
  theme(panel.background = NULL)+
  labs(x="Kd", y= "Frequency")+
  geom_vline(xintercept=2.62)#mean(dib_Kd$Kda)
  
loghl_his <- 
dib_sf |> ggplot(aes(log(vida_media)))+
  #geom_histogram(binwidth = 0.4) +
  geom_density(aes(y = ..density.. * (nrow(dib_hl) * 0.4)), fill="#226E53")+
  #theme_bw()+
  theme(panel.background = NULL)+
  labs(x="Ln(half-life)", y= "Frequency")+
  geom_vline(xintercept=2.515)#mean(log(dib_hl$vida_media))

logkd_his <- 
dib_sf |> ggplot(aes(log(Kda)))+
  #geom_histogram(binwidth = 0.2) +
  geom_density(aes(y = ..density.. * (nrow(dib_Kd) * 0.2)), fill="#FFA700")+
  #theme_bw()+
  theme(panel.background = NULL)+
  labs(x="Ln(Kd)", y= "Frequency")+
  geom_vline(xintercept = 0.818)#mean(log(dib_Kd$Kda))

library(ggpubr)

ggarrange(hl_his,loghl_his,kd_his,logkd_his)

library(ggpmisc)
  
dib_sf |> ggplot(aes(Kda,log(vida_media), col=Zone_4))+
  geom_point()+
  #stat_poly_line() +
  #stat_cor(p.accuracy = 0.001, r.accuracy = 0.05)+
  scale_x_continuous(limits = c(0,6))+
  theme(panel.background = NULL)

dib_sf |> ggplot(aes(log(Kda),log(vida_media)))+
  geom_point() 
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
points(loc.obs, pch=16, col="#FF5300")

## spde #### with prior specification 

spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(10000, 0.01), # P(range < rangemin) = 0.01
  prior.sigma = c(10, 0.01)) # P(sigma > sigmamax) = 0.01

## hyper #### 
hyper <- list(beta = list(prior = 'normal', param = c(0, 1)))

## formula ####

form <- y ~ -1 +
  intercept_Kd +
  intercept_hl + 
  
  SOC_kd+
  #SOC_hl+
  clay_kd+
  land_use_hl+
  
  f(s1, model = spde) + 
  f(s2, model = spde) + 
  f(s12, copy = "s1", fixed = FALSE, hyper = hyper) 

## projector matrices ####


Akd <- inla.spde.make.A(mesh, loc_Kd) 
Ahl <- inla.spde.make.A(mesh, loc_hl) 

#A <- inla.spde.make.A(mesh, loc.obs)

## stacks ####

stackkd <- inla.stack(
  data = list(y = cbind(NA, log(dib_Kd$Kda))),
  A = list(Akd, 1), 
  effects = list(
    list(s1 = 1:spde$n.spde),
    list(intercept_hl  = rep(1, nrow(dib_Kd)),
         SOC_hl=rep(NA,nrow(dib_Kd)),
         SOC_kd=dib_Kd$SOC,
         land_use_hl=rep(NA,nrow(dib_Kd)),#as.character(dib_Kd$adapt_obs),
         clay_kd=dib_Kd$Clay
         )
    ),
  tag = 'kd')

stackhl <- inla.stack(
  data = list(y = cbind(log(dib_hl$vida_media), NA)),
  A = list(Ahl, 1), 
  effects = list(
    list(s1 = 1:spde$n.spde, s12 = 1:spde$n.spde),#s12 = 1:spde$n.spde
    list(intercept_Kd = rep(1, nrow(dib_hl)),
         SOC_hl=dib_hl$SOC,
         SOC_kd=rep(NA,nrow(dib_hl)),
         land_use_hl=as.factor(dib_hl$adapt_obs),
         clay_kd=rep(NA,nrow(dib_hl))#scale(dib_hl$Clay)
        )
    ),
  tag = 'hl') 



stack <- inla.stack(stackhl, stackkd) 


# hyper.eps <- list(hyper = list(prec = list(prior = 'pc.prec',
#                                   param = c(10000, 0.01))))

### empirical data ####
form <- y ~ -1 +
  intercept_hl + 
  intercept_Kd +
  
  SOC_kd+
  #SOC_hl+
  clay_kd+
  land_use_hl+
  
  f(s1, model = spde) + 
  f(s2, model = spde) + 
  f(s12, copy = "s1", fixed = FALSE, hyper = hyper) 



# e.sd <- c(0.75,0.5,1,1,1,0.5,1,1)
# range <- c(70000,30000)
# mvar <- c(0.2,0.8)
# beta <- c(2.5,0.8,0,0,0,0.5,0,0)
# 
# theta.ini <- c(log(1 / e.sd^2),
#                c(log(range),
#                  log(sqrt(mvar))), beta)
# fitting ####

result <- inla(form, rep('gaussian', 2), 
               data = inla.stack.data(stack),
               #lincomb = lincomb.predict,
               #control.family = list(hyper.eps, hyper.eps), 
               control.predictor = list(compute = TRUE,A = inla.stack.A(stack)),
               #control.mode = list(theta = theta.ini, restart = TRUE),
               #control.inla = list(int.strategy = 'eb'),
               control.fixed = list(#prec=list(intercept_Kd= 1/0.75,
                                              #intercept_hl= 1/0.5),
                                    mean=list(intercept_Kd= 2,
                                              intercept_hl= 0.8)), 
               control.compute = list(return.marginals=TRUE,
                                      return.marginals.predictor=TRUE,
                                      hyperpar=TRUE
                                      )
               )


summary(result)

result$summary.fitted.values


sigmas1 <- inla.tmarginal(function(x) 1/sqrt(exp(x)),result$internal.marginals.hyperpar[[3]])
sigmaepsilon <- inla.tmarginal(function(x) 1/sqrt(exp(x)),result$internal.marginals.hyperpar[[1]])
restab=sapply(result$marginals.fixed, function(x) inla.zmarginal(x,silent=TRUE))
restab=cbind(restab, "alpha" = inla.zmarginal(sigmas1,silent=TRUE))
restab=cbind(restab, "epsilon" = inla.zmarginal(sigmaepsilon,silent=TRUE))
data.frame(restab)

ddf <- data.frame(rbind(sigmas1, sigmaepsilon),
                  errterm = gl(2, nrow(sigmaepsilon), labels = c("alpha", "epsilon")))
ggplot(ddf, aes(x, y, linetype = errterm)) + geom_line() + xlab(" ") +
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


book.plot.field <-
  function(field,
           mesh,
           projector,
           xlim,
           ylim,
           dims = c(300, 300),
           poly,
           asp = 1,
           axes = FALSE,
           xlab = '',
           ylab = '',
           title = '',
           col = book.color.c(),
           polyBoundary = NULL,
           return_sf = FALSE,
           ...) {
    # browser()
    ## you can supply field as a matrix vector or like a named list with 'x', 'y' and 'z' as for image
    ## when field is a vector, it will project it using projector, assuming projector will create a matrix
    ## when mesh is supplied and projector not, projector will be created and used to project field
    if (missing(mesh)) {
      if (missing(projector)) {
        if (missing(xlim) | missing(ylim)) {
          fields::image.plot(
            field,
            asp = asp,
            axes = axes,
            xlab = xlab,
            ylab = ylab,
            col = col,
            ...
          )
        } else {
          fields::image.plot(
            field,
            xlim = xlim,
            ylim = ylim,
            asp = asp,
            axes = axes,
            xlab = xlab,
            ylab = ylab,
            col = col,
            ...
          )
        }
      } else {
        if (missing(xlim))
          xlim <- range(projector$x)
        if (missing(ylim))
          ylim <- range(projector$y)
        field.proj <- INLA::inla.mesh.project(projector, field)
        
        if (is(field.proj, 'matrix')) {
          # Convert matrix + projector to raster.
          raster <-
            raster::raster(t(field.proj)[nrow(field.proj):1,])
          raster::extent(raster) <-
            c(range(projector$x), range(projector$y))
          
          if (!is.null(polyBoundary)) {
            r2 <- raster::crop(raster, raster::extent(polyBoundary))
            raster <- raster::mask(r2, polyBoundary)
          }
        }
        
        
        x <- y <- NULL
        # browser()
        raster.df <- methods::as(raster, "SpatialPixelsDataFrame")
        sf_df <- st_as_sf(methods::as(raster, "SpatialPixelsDataFrame"))
        st_crs(sf_df) <- st_crs(polyBoundary)
        raster.df <- as.data.frame(raster.df)
        raster.df <-
          tidyr::gather(raster.df, key = "raster_name", value = "z", -x, -y)
        
        names(raster.df) <- c('long', 'lat', 'raster_name', 'value')
        
        
        plot <- ggplot2::ggplot() +
          ggplot2::geom_raster(data = raster.df,
                               aes_string('long', 'lat', fill = 'value')) +
          ggplot2::labs(x = xlab, y = ylab, title  = title) +
          ggplot2::theme_void()
        if (return_sf) {
          return(list(plot = plot,
                      sf_obj = sf_df))
        } else {
          return(plot)
        }
        
      }
    } else {
      if (missing(xlim))
        xlim <- range(mesh$loc[, 1])
      if (missing(ylim))
        ylim <- range(mesh$loc[, 2])
      projector <- INLA::inla.mesh.projector(mesh,
                                             xlim = xlim,
                                             ylim = ylim,
                                             dims = dims)
      field.proj <- INLA::inla.mesh.project(projector, field)
      
      
      fields::image.plot(
        x = projector$x,
        y = projector$y,
        z = field.proj,
        asp = asp,
        axes = axes,
        xlab = xlab,
        ylab = ylab,
        col = col,
        ...
      )
    }
    if (!missing(poly))
      plot(poly, add = TRUE, col = 'grey')
  }

# Posterior mean of 'm'
par(mfrow = n2mfrow(length(names(result$summary.ran)), asp = 1))
# Posterior mean of 'm'
for (i in names(result$summary.ran)) {
  book.plot.field(field = result$summary.random[[i]]$mean, projector = prj, main = i) 
}

maps
s1 <- book.plot.field(field = result$summary.random[["s1"]]$mean, 
                projector = prj, 
                polyBoundary = limits_sf,
                return_sf = T,
                main = i) 


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
