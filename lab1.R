#
# print('Random Proceses Lab1 Nedilko Daryna', quote=FALSE)
#
# # generating random numbers and building a plot
n <- 1000000
# dots <- runif(n)
# # The size of graph field
# par(mar=c(1, 1, 1, 1))
# plot(dots, pch = 1, cex = 0.01 )
# #
# #
# expected_value <- sum(dots/n)
# cat('Expected value ', expected_value,"\n")
# variance <- sum(dots^2)/n - expected_value^2
# cat('Variance: ',variance,"\n")
# h <-1/(1+3.322*log10(n))
# cat("Interval length: ",h,"\n")
# nint <- round(1/h)
# cat('Number of intervals: ', nint,"\n")
#
# # Creating intervals, calculaving frequensy
# bounders <- c(seq(0, 1, h),1)
# # print(range(dots))
# l <- as.numeric(length(bounders))-1
# freq <- c()
# for (i in 1:l){
#   suB <- 0
#   for (k in dots){
#     if(k>=bounders[i] && k<bounders[i+1]){
#       suB <- suB+1}}
#   freq <- c(freq,suB)
# }
#
# intervals <-c()
# for (i in 1:l){
#   temp <- paste('[',as.character(round(bounders[i],digits=5)),',',as.character(round(bounders[i+1],digits=5)),')')
#   intervals <- c(intervals,temp)
# }
# frec_comp <- freq/n
# # table.Inter <- intervals
# # table.freq <- freq
# # table.frec_comp <- frec_comp
# table <- cbind(intervals, "|",round(freq,digits=5),"|",round(frec_comp,digits=5))
# print(table, quote=FALSE )
#
# # Building a histogram
# library(ggplot2)
# fhist <- cut(dots,bounders, right=FALSE)
# # barplot(table(fhist),main="Histogram for uniform distribution of dots", col = "white", border = "steelblue" )
#
# # hyposesys testing for unified distribution
# # будемо використовувати критерій хі квадрат
# # спочатку обрахуємо середньоарифметичне кожного з інтервалів
# evgper <- c()
# for (i in 1:l){
#   avg <- (bounders[i]+bounders[i+1])/2
#   evgper <- c(evgper,avg)
# }
# print(evgper)
# # виводимо таблицю середнього значення по інтервалу та частоту попадання у цей інтервал.
# table1 <- print(cbind("|",evgper,"|",freq), quote = FALSE)
# #Обчислюємо середьоквадратичне відхилення, вибіркове середнє нам вже відомо
# standart_deviation <- sqrt(sum(dots-expected_value)^2/(n-1))
# cat("Expected value: ", expected_value,"\n")
# cat("Standart deviation",standart_deviation,"\n")
#
# # Оцінка параметрів a,b - кінців інтервалу
# adot <- expected_value - sqrt(3)*standart_deviation
# bdot <- expected_value + sqrt(3)*standart_deviation
# # Знаходимо щільність гіпотетичного рівномірного розподілу
# fx <- 1/(bdot-adot)
# cat("Hypotetic frequensy of unified distribution", fx, "\n")
#
# # Theoretical friquensie
# n1 <- n*fx*(evgper[1]-adot)
# nl <- n*fx*(bdot-evgper[l])
# nth <- c(n1)
# for (i in 2:(l-1)){
#   ni <- n*(evgper[i]-evgper[i+1])*fx
#   nth <-c(nth,ni)
# }
# nth <- c(nth,nl)
# print(nth)
# forHiSpos <- c()
# for (i in 1:l){
#   calc <- ((freq[i]-nth[i])^2)/nth[i]
#   forHiSpos <- c(forHiSpos,calc)
# }
#
# # Використовуємо критерій Пірсона при k=l-3
# # зводимо усі дані до таблиці
# table2 <- cbind(intervals,"|",round(evgper,digits = 5),"|",round(nth,digits=6),"|",round(forHiSpos,digits=6))
# print(table2)
# # знаходимо хі квадрат
# HiSQRTspos <- sum(forHiSpos)
# # при рівні значущості 0.95 - введено вручну, тому для n відміннне від мільйона треба зазначати інше критичне значення
# HiCrit <- 7.01
# # Перевірка чи гіпотеза правильна
# if (HiSQRTspos<HiCrit){
#   cat("Hyposesys about unified distribution IS TRUE")
# }else{cat("That`s not unified distribution")}
# #
# #
# 1.4
# Дослідити розподіл випадкової величини
# maxes <- c()
# m <- 1000
# for (i in 1:m){
#   tempmax <- as.numeric(max(runif(n)))
#   max <- c(max,tempmax)
# }
# hist(maxes)


# #2.0
# xi42 <- c(1,10,15,23,29,38,42)
# pi42 <- c(0.02, 0.05,0.1, 0.28, 0.23, 0.22,0.1)
#
# table21 <- print(cbind(xi42,"|", pi42), quote=FALSE)
#
# # generating set
# nintervals <- length(pi42)
# discrete_dist <- sample(x=xi42,n,replace = T, prob = pi42)
# hist(discrete_dist)
# varDiscDest <- var(discrete_dist)
# MathExpectation <- mean(discrete_dist)
# cat("Variance for this disctere distribution: ",varDiscDest,"\nExpected value: ", MathExpectation)
# forPlot <- cut(discrete_dist, c(xi42, 100), right = FALSE)
# barplot(table(forPlot),main="Histogram for discrete", col = "white", border = "steelblue" )
# print(table(forPlot))
# freq2 <- c()
# for (i in xi42){
#   suB2 <-0
#   for (k in discrete_dist){
#     if (i==k){suB2 <- suB2+1 }
#   }
#   freq2 <- c(freq2,suB2)
# }
# print(freq2)
#
# #
# stand_dev42 <- sqrt(varDiscDest)
# cat("Standart deviation: ", stand_dev42, "\nExpected value: ", MathExpectation)
# h <- (max(xi42)-min(xi42))/length(xi42)
#
#
# # gauss <- function(u){
# #   g <- (2.718281828459045^(-u^2/2))/sqrt(2*pi)
# #   return (g)
# # }
# k <- print(sum(((freq2-n*pi42)^2)/n*pi42))
# # ui
#
# 4
# a1<-0
# b1<-1
# I1 <- print((b1-a1)*sum(sapply(dots,function (x) x^7+x^5+x^3 ))/n)
#
# a2<-0
# b2<-pi
# newdots <- runif(n,a2,b2)
# I2 <- print((b2-a2)*sum(sapply(newdots,function (x) 2*sin(3*x) ))/n)
#
# a3<-0
# b3<-1000000
# newwdots <- c(0:1000000)
# I3 <- print((sum(sapply(newwdots, function (x) 1/((x+1)*sqrt(x)) ))*(b3-a3)/n))
#

# bet <- c()
# dots <- rnorm(n)
# # hist(fb)
# # sapply
#
# bet<- sapply(c(1:n), function(x) +(dots[x+1]<dots[x]))
#
# hist(bet, breaks = c(0,0.4,0.6,1) )

# function(a,b,f, dots){
#   print(something)
#   INT <- (b-a)*sum(f(dots))/length(dots)
#   return (INT)
# }

# hist(bounders,freq)
# for (i in 1:)

# hist(dots, breaks = bounders, prob=TRUE)
# ranges <- c(cut(dots,bounders, right=FALSE))
# print(ranges)

# freq <- print(table(dots.cut))
# tab <- print(matrix(ranges))
# colnames(tab) <- c("Interval", "Frequensy","Comparable frequensy")
# tab <- as.table(tab)

# table(dots.cut),table(dots.cut)/n



# Monte-Carlo method
# install.packages('plotrix')
# library(plotrix)
# n<-15
# draw.circle(0,0,2)

# plot(0,0,type = "n",xlim = c(-2*n,2*n), ylim = c(-2*n,2*n))
# xs <- sample(-n:n, size = n)
# ys <- sample(-n:n, size = n)
# rs <- sample(0:n)
# cat('Table of circles')
# for_t <- print(cbind('|',xs,'|',ys,'|',rs,'|'), quote=FALSE)
# for(i in 1:n){
#   draw.circle(xs[i],ys[i],rs[i])
# }
windows()
library(rgl)
n<-as.numeric(readline("Enter the amount of spheres"))
g<-n/2
# plot(0,0,type = "n",xlim = c(-2*n,2*n), ylim = c(-2*n,2*n))
xs <- sample(-g:g,size=n)
ys <- sample(-g:g,size=n)
zs <- sample(-g:g,size=n)
rs <- sample(g:n, size=n, replace = TRUE)
cat('Table of spheres\n')
for_t <- print(cbind('|',xs,'|',ys,'|',zs,'|',rs,'|'), quote=FALSE)
cords <- data.frame(xs,ys,zs)
open3d()
material3d(1,alpha=0.2)
spheres3d(xs,ys,zs,rs,color = rainbow(n))
# for(i in 1:n){
#   draw.circle(xs[i],ys[i],rs[i])
# }
#
# library(ggplot2)
# # ggplot(data = dots_for_6, aes(x = weight, y = hindfoot_length)) + geom_point()
# # newn<-1000000
# # dots_for_6 <- matrix(runif(2*newn,-2*n,2*n),nrow=newn)
# # plot(dots_for_6[,1],dots_for_6[,2],pch = 1, cex = 0.1)
check<- function (x){
  for (i in x){
    if (i == FALSE){
      return (FALSE)}
  }
  return (TRUE)
}
#
#
# # func_to_c <- function(xs,ys,i){
# #   temp<-sapply(1:n,function(k) sqrt((i[1]-xs[k])^2 + (i[2]-ys[k])^2)<=rs[k])
# #   if (check(temp)==TRUE){
# #     return (c(i[1],i[2]))
# #   }
# # }
#
# do_intersec <- check(circles.intersect(cords,rs))
# if(check(circles.intersect(cords,rs))==FALSE){cat("Circles don`t intersect")}else{cat("Calculating intersection area")}
#
#
# # inersection wiil be in smallest circle, so w
# #
# # e can fill only sq around this circle
v<-min(rs)
find_small <- which(rs==v)[1]
xss<- xs[find_small]
yss <- ys[find_small]
zss <- zs[find_small]
# #
# # print(cords[1][5])
# # calculating num of dots in intersection of circles
#
#
#

newn<-1000000
dots_for_6 <- matrix(c(runif(newn,(xss-v),(xss+v)),runif(newn,(yss-v),(yss+v)),runif(newn,(zss-v),(zss+v))), ncol=3)
kz <- sapply(1:newn, function(i) if(check(sapply(1:n,
                                        function(k) sqrt((dots_for_6[i,1]-xs[k])^2 + (dots_for_6[i,2]-ys[k])^2+(dots_for_6[i,3]-zs[k])^2)<rs[k]))){dots_for_6[i,1:3]})

full<- print(is.null(unlist(Filter(Negate(is.null),kz))))
if (full){
  cat("This spheres have no intersection")
}else{
  dots_for_xyz<-t(matrix(unlist(Filter(Negate(is.null),kz)),nrow=3))
  spheres3d(xs,ys,zs,rs,color = rainbow(n))
  points3d(dots_for_xyz[,1],dots_for_xyz[,2],dots_for_xyz[,3], color="black", alpha=1)

  volume<-(((2*v)^3)*length(dots_for_xyz))/newn
  cat("The volume of intersection is: ",volume)

}

# circles.plot(cords,rs)
# points(dots_for_6,col="green")
# points(dots_for_xy)
# points(dots_for_6,pch = 1, cex = 0.01,col="green")



# print(dots_new_xy)
# dots_matrix <- print(matrix(dots_new_xy,ncol=2))
# plot(, asp = 1, xlim = c(-2*n, 2*n),ylim = c(-2*n, 2*n))
# print(dots_new_y)
# # print(dots_new_x)
# l<-length(dots_for_xy)/2
# all_graph<- rbind(cords,dots_for_xy)
#
#
# circles.plot(dots_for_xy,)
# draw.circle(xs,ys,rs)
# data_for_pic <- print(data.frame(dots_new_x,dots_new_y))
# draw.circle(dots_new_x,dots_new_y,0)



















# check_in_intersec <-c()

#
