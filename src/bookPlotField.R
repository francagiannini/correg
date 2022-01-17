book.plot.field <- function(field, mesh, projector, xlim, ylim, 
                            dims=c(300,300), poly, asp = 1, 
                            axes = FALSE, xlab = '', ylab = '', 
                            col = book.color.c(), ...){
  ## you can supply field as a matrix vector or like a named list with 'x', 'y' and 'z' as for image
  ## when field is a vector, it will project it using projector, assuming projector will create a matrix 
  ## when mesh is supplied and projector not, projector will be created and used to project field
  if (missing(mesh)) {
    if (missing(projector)) {
      if (missing(xlim) | missing(ylim)) {
        fields::image.plot(field, asp = asp, axes = axes, 
                           xlab = xlab, ylab = ylab, col = col, ...)
      } else {
        fields::image.plot(field, xlim = xlim, ylim = ylim, asp = asp, 
                           axes = axes, xlab = xlab, ylab = ylab, col = col, ...)
      }
    } else {
      if (missing(xlim)) xlim <- range(projector$x)
      if (missing(ylim)) ylim <- range(projector$y)
      field.proj <- inla.mesh.project(projector, field)
      fields::image.plot(x = projector$x, y = projector$y, z = field.proj, 
                         asp=asp, axes=axes, xlab = xlab, ylab = ylab, 
                         col=col, xlim=xlim, ylim=ylim, ...)
    }
  } else {
    if (missing(xlim)) xlim <- range(mesh$loc[,1])
    if (missing(ylim)) ylim <- range(mesh$loc[,2])
    projector <- inla.mesh.projector(mesh, xlim = xlim,
                                     ylim = ylim, dims=dims)
    field.proj <- inla.mesh.project(projector, field)
    fields::image.plot(x = projector$x, y = projector$y, z = field.proj, 
                       asp=asp, axes=axes, xlab = xlab, ylab = ylab, col=col, ...)
  }
  if (!missing(poly)) 
    plot(poly, add = TRUE, col = 'grey')
}


book.color.c = function(n = 201) {
  return(viridis::viridis(n))
}
