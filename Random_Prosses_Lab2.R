# windows()

pdf(file = "C:/Users/Meduzka/Desktop/My Plot.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches


n <- 100
events <- 1:n
rand <- runif(n)
l <- as.numeric(unlist(readline('Enter lambda(intensivity)')))
Sk <- -log(1-rand)/l
time <- sapply(events,function(k) sum(Sk[1:k]))
# windows()
# par(mfcol=c(5,3),mai=c(0.5,0.5,0.5,0))
plot(events,time, pch=20)
# '''-----'''
T <- n/l
Ts <- sapply(events, function(k) k/l)
library(ggplot2)
df_1 <- data.frame(events,Ts)
ggplot(df_1, aes(x = events, y = Ts)) +
       geom_bar(stat="identity")
# '''-----'''
df_2 <- data.frame(events,Sk)
ggplot(df_2, aes(x = events, y = Sk)) +
       geom_bar(stat="identity")
# '''-----'''
Pn <- exp(-l*T)*((l*T)^n)/factorial(n)
cat('Frequensy of n events aproximately equal ',Pn)

repetition <- 1:(n*11)
ffreq <- sapply(repetition, function(i) sum(sapply(runif(n),function(k) -log(1-k)/l)))*l
frequency <- sapply(ffreq, function(e) exp(-e)*(e^n)/factorial(n))

df_3 <- data.frame(repetition,frequency)
ggplot(df_3, aes(x = repetition, y = frequency)) + ylim(0,max(frequency)+0.001) + geom_bar(stat="identity") +
  geom_hline(yintercept=Pn, linetype="dashed", color = "red")

custom <- function(x){
  if (round(x,3)==round(Pn,3)){
    return (1)
  }else{
    return (0)
  }

}


events_n <- sum(sapply(frequency, function(i) custom(i)))
cat('\nFrequency of n events amount 1100 Poisson prosses models is ', events_n)

# dev.off()

# ---------------вінерівський випадковий процес


# library('reshape2')

W1 <- function (M,t){
  mu <- runif(M+1)
  t <- runif(1)
  value <- t*mu[1] + sqrt(2) * sum(sapply(1:M,function(i) sin(i*pi*t)*mu[i+1]/(pi*i)))
  return (value)
}

W2 <- function(M,t){
  mu0 <- runif(1)
  mu1 <- runif(M)
  mu2 <- runif(M)
  t <- runif(1)
  value <- sapply(1:M, function(i) t*mu0 + sqrt(2) * sum(mu1[i]*sin(2*pi*i*t)/(2*pi*i)+mu2[i]*(1-cos(2*pi*i*t))/(2*pi*i)))
  return (value)
}

Time <- 100
M1 <- 1000
M2 <- 10000
M3 <- 100000
t <- 1:Time
m1_w1 <- sapply(t,function(i) W1(M1,i))
m1_w2 <-sapply(t,function(i) W2(M1,i))
m2_w1 <- sapply(t,function(i) W1(M2,i))
m2_w2 <-sapply(t,function(i) W2(M2,i))
m3_w1 <-sapply(t,function(i) W1(M3,i))
m3_w2 <-sapply(t,function(i) W2(M3,i))
# df_m1 <- data.frame(t, m1_w1, m1_w2)
# df_m2 <- data.frame(t, m2_w1, m2_w2)
# df_m3 <- data.frame(t, m3_w1, m3_w2)
df_m1_w1 <- data.frame(t, m1_w1)
colnam1 <- c('Time', 'W1')
df_m1_w2 <- data.frame(t, m1_w2)
colnam2 <- c('Time', 'W2')
colnames(df_m1_w1)= colnam1
colnames(df_m1_w2)= colnam2
print('---------------------------------------------------------------------------')
# df_m1.m <- melt(df_m1, id.vars = 't')
# df_m2.m <- melt(df_m2, id.vars = 't')
# df_m3.m <- melt(df_m3, id.vars = 't')
df_m2_w1 <- data.frame(t, m2_w1)
df_m2_w2 <- data.frame(t, m2_w2)
colnames(df_m2_w1)= colnam1
colnames(df_m2_w2)= colnam2

df_m3_w1 <- data.frame(t, m3_w1)
df_m3_w2 <- data.frame(t, m3_w2)
colnames(df_m3_w1)= colnam1
colnames(df_m3_w2)= colnam2

ggplot() + geom_line(data=df_m1_w1, aes(x=Time, y=W1),color='red',shape=18) +geom_line(data=df_m1_w2, aes(x=Time, y=W2),color='blue',shape=20)+ ylab('W2 and W1')

ggplot() + geom_line(data=df_m2_w1, aes(x=Time, y=W1),color='red',shape=18) +geom_line(data=df_m2_w2, aes(x=Time, y=W2),color='blue',shape=20)+ ylab('W2 and W1')

ggplot() + geom_line(data=df_m3_w1, aes(x=Time, y=W1),color='red',shape=18) +geom_line(data=df_m3_w2, aes(x=Time, y=W2),color='blue',shape=20)+ ylab('W2 and W1')

# ggplot(df_m2.m, aes(t,value, color = variable)) + geom_point(shape=18) + scale_color_manual(values=c('red','blue'))
# ggplot(df_m3.m, aes(t,value, color = variable)) + geom_point(shape=18) + scale_color_manual(values=c('red','blue'))
mean_W1 <- c(mean(m1_w1), mean(m2_w1),mean(m2_w1))
mean_W2 <- c(mean(m1_w2), mean(m2_w2),mean(m2_w2))
variance_W1 <- c(var(m1_w1), var(m2_w1),var(m2_w1))
variance_W2 <- c(var(m1_w2), var(m2_w2),var(m2_w2))
MS <- as.integer(c(M1,M2,M3))
mean_M <- table(c(MS,mean_W1,mean_W2))
# var_M <- matrix(c(MS,variance_W1,variance_W2),3,3)
# cat("\n Table of mean values\n M  W1  W2\n")
for (i in 1:3){
  cat("\n For function W1 where M= ",MS[i],"\nmean of distibution is ", mean_W1[i],"\nvariance is ", variance_W1[i])
  cat("\n For function W2 where M= ",MS[i],"\nmean of distibution is ", mean_W2[i],"\nvariance is ", variance_W2[i])
}
