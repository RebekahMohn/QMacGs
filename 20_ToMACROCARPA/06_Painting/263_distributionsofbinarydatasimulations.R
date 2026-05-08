library(ggplot2)
library(ggpubr)

ex1<-rbinom(n=375,size=5,prob=0.5)/5

ex2<-rbinom(n=375,size=10,prob=0.5)/10

ex3<-rbinom(n=375,size=20,prob=0.5)/20

ex4<-rbinom(n=375,size=40,prob=0.5)/40

ex5<-rbinom(n=375,size=100,prob=0.5)/100

ex6<-rbinom(n=375,size=400,prob=0.5)/400

ggarrange(
  ggplot()+
    geom_histogram(aes(x=ex1))+
    geom_vline(xintercept = qbinom(8.699813e-09,size=5,prob=.5)/5)+
    geom_vline(xintercept = qbinom(1-8.699813e-09,size=5,prob=.5)/5)+
    xlim(0,1),
  ggplot()+
    geom_histogram(aes(x=ex2))+
    geom_vline(xintercept = qbinom(8.699813e-09,size=10,prob=.5)/10)+
    geom_vline(xintercept = qbinom(1-8.699813e-09,size=10,prob=.5)/10)+
    xlim(0,1),
  ggplot()+
    geom_histogram(aes(x=ex3))+
    geom_vline(xintercept = qbinom(8.699813e-09,size=20,prob=.5)/20)+
    geom_vline(xintercept = qbinom(1-8.699813e-09,size=20,prob=.5)/20)+
    xlim(0,1),
  ggplot()+
    geom_histogram(aes(x=ex4))+
    geom_vline(xintercept = qbinom(8.699813e-09,size=40,prob=.5)/40)+
    geom_vline(xintercept = qbinom(1-8.699813e-09,size=40,prob=.5)/40)+
    xlim(0,1),
  ggplot()+
    geom_histogram(aes(x=ex5))+
    geom_vline(xintercept = qbinom(8.699813e-09,size=100,prob=.5)/100)+
    geom_vline(xintercept = qbinom(1-8.699813e-09,size=100,prob=.5)/100)+
    xlim(0,1),
  ggplot()+
    geom_histogram(aes(x=ex6))+
    geom_vline(xintercept = qbinom(8.699813e-09,size=400,prob=.5)/400)+
    geom_vline(xintercept = qbinom(1-8.699813e-09,size=400,prob=.5)/400)+
    xlim(0,1)
  
)
