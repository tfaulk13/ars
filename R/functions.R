


library(pracma)

params.r = 2
params.m = 10
params.mu = 0
params.sig2 = 1
xmin = -4
xmax = 4

h = function(y){
  v = params.r*y - params.m * log(1+exp(y)) - (y-params.mu)^2/(2*params.sig2) # plus normalizing const
  return(v)
}

## derivative of h (replace this) ---------------------------------------------
#
#
#
#
dh = function(y){
  params.r - params.m * exp(y) / (1 + exp(y)) - (y-params.mu)/params.sig2
}
x = seq(-6, 6, length.out = 100)
plot(x, h(x), type="l")
```
```{r}


## z_k ------------------------------------------------------------------------
#
#
#
#
find_intercept = function(x_pts){
  # x_pts: points selected initially to find the tangent line, abscissae, T_k
  # z: intercept x values
  x0 = head(x_pts, n=-1)
  x1 = tail(x_pts, n=-1)
  z = (h(x1) - h(x0) - x1*dh(x1) + x0*dh(x0)) / (dh(x0) - dh(x1))
  return(z)
}
# g1 = function(x) x^(5-1)*exp(-2*x)
# g2 = function(x) exp(-x^2/2)
# h2 = function(x) log(exp(-x^2/2))
# fderiv(f, 4, n = 1, h = 0)
a = c(-3, -0.2, 2)
find_intercept(a)



# u_k, the upper bound --------------------------------------------------------
#
#
#
#
find_u = function(x, x_pts) {
  # x: points randomly selected to plot the bounds
  # x_pts: points selected initially to find the tangent line, abscissae, T_k
  # u: the corresponding y values for x
  # norm_const: area under the curve

  z_pts = find_intercept(x_pts) #length(z_pts) = length(x_pts)-1
  all_pts = c(xmin, z_pts, xmax) # if connect all_pts, can form the upper hull
  u = c() # u is value of u_k function, the y value at each of the all_pts
  norm_const = c()
  # find which piece is each x in
  x_idx = findInterval(x, c(xmin, z_pts, xmax))
  # loop through all pieces, for each piece, calc the corresponding y for each x
  for(i in 1:(length(all_pts)-1)){
    xp = x[x_idx == i]
    uu = h(x_pts[i]) + (xp - x_pts[i])*dh(x_pts[i])
    u[x_idx == i] = uu
    # find the area under each piece
    norm_const[i] = integrate(function(z) h(x_pts[i]) + (z - x_pts[i])*dh(x_pts[i]), lower=all_pts[i], upper = all_pts[i+1])[[1]]
  }
  # u[1] = h(x_pts[1]) + (xmin-x_pts[1])*dh(x_pts[1])
  # # loop through all the other all_pts, start from the second number
  # for (i in 2:(length(all_pts))){
  #   u[i] = h(x_pts[i-1]) + (all_pts[i]-x_pts[i-1])*dh(x_pts[i-1])
  #   norm_const[i-1] = integrate(function(z) exp(h(x_pts[i-1]) + (z-x_pts[i-1])*dh(x_pts[i-1])), lower=all_pts[i-1], upper = all_pts[i])[[1]]
  # }
  #return(u)
  return(list(u=u, norm_const=norm_const))
}



# s_k -------------------------------------------------------------------------
#
#
#
#
find_s = function(x, x_pts){
  # x: points randomly selected to plot the bounds
  # x_pts: points selected initially to find the tangent line, abscissae, T_k

  s = exp(find_u(x, x_pts)$u)/sum(find_u(x, x_pts)$norm_const)
  return(s)
}
find_s(c(-1, 0, 1), a)


# l_k, lower bound ------------------------------------------------------------
#
#
#
#
find_l = function(x, x_pts){
  # x: points randomly selected to plot the bounds
  # x_pts: points selected initially to find the tangent line, abscissae, T_k
  # l: the corresponding y values for x

  x_idx = findInterval(x, x_pts)
  l = c()
  for(i in 1:(length(x_pts)-1)){
    xp = x[x_idx == i]
    ll = ((x_pts[i+1] - xp)*h(x_pts[i]) + (xp-x_pts[i])*h(x_pts[i+1])) / (x_pts[i+1]-x_pts[i])
    l[x_idx == i] = ll
  }
  return(l)
}



# test: plot all functions created
a = c(-3, -0.2, 2)
x = seq(-4, 3.9, length.out = 500)
plot(x, h(x), type="l", col="red", ylim=c(-50, 0))
lines(x, find_l(x, a))
lines(x, find_u(x, a)$u)

