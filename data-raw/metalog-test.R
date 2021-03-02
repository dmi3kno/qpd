x <- rmetalog::fishSize$FishSize
step_len=0.01

l <- length(x)


make_sym_grid <- function(l){
  # from Gilchrist(2000), p5
  (seq(l)-0.5)/l  
}


make_finetailed_grid <- function(x, step_len = 1/length(x)){
  tailstep <- step_len / 10
  
  c(seq(tailstep, by=tailstep, length.out = 9),
    seq(step_len, 1-step_len, step_len),
    seq(from=(1 - step_len + tailstep), by=tailstep,  length.out = 9))
}


x_new <- stats::quantile(x, probs = y_new)
length(x_new)

plot(y~x, type="l")
lines(y_new~x_new, col=2)


y <- seq(step_len, (1 - step_len), step_len)
tailstep <- (step_len / 10)

seq(tailstep, (min(y) - tailstep), tailstep)
