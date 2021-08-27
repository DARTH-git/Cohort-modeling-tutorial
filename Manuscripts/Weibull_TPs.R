hweibull(3, 2, 1)
dweibull(3, 2, 1) / pweibull(3, 2, 1, lower.tail = FALSE)
Hweibull(3, 2, 1)
-pweibull(3, 2, 1, log.p = TRUE, lower.tail = FALSE)
v_time <- 1:75
scale = 0.08 #lambda
shape = 1.20
H <- (v_time*scale)^shape
plot(H)
H_u <- ((v_time-1)*scale)^shape

h <- shape*scale*(v_time*scale)^(shape-1)
plot(h)
points(1-exp(-h), col = "blue")

tp <- 1-exp(H_u - H)
tp
points(tp, col = "red")

# tp_old <- scale*shape*v_time^(shape-1)
points(tp_old, col = "green")

tp_new <- 1-exp(((v_time-1)*scale)^shape - (v_time*scale)^shape)
