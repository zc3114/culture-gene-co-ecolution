library("reshape2")
library("ggplot2")

f <- function(b0=0,b1=0.5,b2=0.5,b3=1,d=0.5,s=0,compen=3.68,pop=80000,alpha=0,f=0,gen=1440){
  A2=1/2/pop  #starting frequencies of mutation allele 
  A1=1-A2
  A1A1=A1*A1
  A1A2=A1*A2*2
  A2A2=A2*A2
  #allele frequencies today
  a2=0.916
  a1=1-a2
  
  h=1-1/pop
  
  ng <- function(x1,x2,x3,x4,x5,x6,alpha,d,b0,b1,b2,b3,s,compen,f){ #s=true selection coefficient comp=compensation factor
    s<-s/compen
    z<-1+alpha*d*x3+alpha*d*x4+alpha*x5+alpha*x6; #normalising factor
    a<-(1+alpha)
    da<-(1+d*alpha)
    #next gen pheno-genotype, without selection
    A1A1h<-b3*(x1*x1+x1*x3+1/4*x3*x3)+b2*(x1*x2+0.5*x2*x3+0.5*x1*x4+0.25*x3*x4)+b1/z*(x2*x1+0.5*x2*x3*da+0.5*x1*x4+0.25*da*x3*x4)+b0/z*(x2*x2+0.5*x2*x4*da+0.5*x4*x2+0.25*x4*x4*da)
    A1A1H<-(1-b3)*(x1*x1+x1*x3+1/4*x3*x3)+(1-b2)*(x1*x2+0.5*x2*x3+0.5*x1*x4+0.25*x3*x4)+(1-b1)/z*(x2*x1+0.5*x2*x3*da+0.5*x1*x4+0.25*da*x3*x4)+(1-b0)/z*(x2*x2+0.5*x2*x4*da+0.5*x4*x2+0.25*x4*x4*da)
    A1A2h<-b3*(2*x1*x5+x1*x3+x3*x5+0.5*x3*x3)+b2*(x1*x6+x2*x5+0.5*(x1*x4+x2*x3+x3*x4+x3*x6+x4*x5))+b1/z*(x2*x5*a+x1*x6+0.5*(x2*x3*da+x1*x4+x3*x4*da+x4*x5*a+x3*x6*da))+b0/z*(x2*x6*a+x2*x6+0.5*(x2*x4*da+x2*x4+x4*x4*da+x4*x6*a+x4*x6*da))
    A1A2H<-(1-b3)*(2*x1*x5+x1*x3+x3*x5+.5*x3*x3)+(1-b2)*(x1*x6+x2*x5+0.5*(x1*x4+x2*x3+x3*x4+x3*x6+x4*x5))+(1-b1)/z*(x2*x5*a+x1*x6+0.5*(x2*x3*da+x1*x4+x3*x4*da+x4*x5*a+x3*x6*da))+(1-b0)/z*(x2*x6*a+x2*x6+0.5*(x2*x4*da+x2*x4+x4*x4*da+x4*x6*da+x4*x6*a))
    A2A2h<-b3*(x5*x5+x5*x3+1/4*x3*x3)+b2*(x5*x6+0.5*x6*x3+0.5*x5*x4+0.25*x3*x4)+b1/z*(x5*x6*a+0.5*x4*x5*a+0.5*x3*x6*da+0.25*da*x3*x4)+b0/z*(x6*x6*a+0.5*x6*x4*a+0.5*x6*x4*da+0.25*x4*x4*da)
    A2A2H<-(1-b3)*(x5*x5+x5*x3+1/4*x3*x3)+(1-b2)*(x5*x6+0.5*x6*x3+0.5*x5*x4+0.25*x3*x4)+(1-b1)/z*(x5*x6*a+0.5*x4*x5*a+0.5*x3*x6*da+0.25*da*x3*x4)+(1-b0)/z*(x6*x6*a+0.5*x6*x4*a+0.5*x6*x4*da+0.25*x4*x4*da)  
    #adding selection coefficient
    nor<-((1+s)*(A2A2h+A2A2H)+(1+d*s)*(A1A2h+A1A2H)+(A1A1h+A1A1H))#normalising
    A1A1hx<-A1A1h/nor
    A1A1Hx<-A1A1H/nor
    A1A2hx<-A1A2h*(1+d*s)/nor
    A1A2Hx<-A1A2H*(1+d*s)/nor
    A2A2hx<-A2A2h*(1+s)/nor
    A2A2Hx<-A2A2H*(1+s)/nor
    
    #adding oblique transmission
    A1A1H <- A1A1Hx+f*A1A1hx*(A1A1Hx+A1A2Hx+A2A2Hx)
    A1A2H <- A1A2Hx+f*A1A2hx*(A1A1Hx+A1A2Hx+A2A2Hx)
    A2A2H <- A2A2Hx+f*A2A2hx*(A1A1Hx+A1A2Hx+A2A2Hx)
    A1A1h <- A1A1hx*(1-f*(A1A1Hx+A1A2Hx+A2A2Hx))
    A1A2h <- A1A2hx*(1-f*(A1A1Hx+A1A2Hx+A2A2Hx))
    A2A2h <- A2A2hx*(1-f*(A1A1Hx+A1A2Hx+A2A2Hx))
    
    return(list("A1A1h"=A1A1h,"A1A1H"=A1A1H,"A1A2h"=A1A2h,"A1A2H"=A1A2H,"A2A2h"=A2A2h,"A2A2H"=A2A2H))
  }
  
  diff = 1
  
  A1A1h=A1A1*h;
  A1A1H=A1A1*(1-h);
  A1A2h=A1A2*h;
  A1A2H=A1A2*(1-h);
  A2A2h=A2A2*h;
  A2A2H=A2A2*(1-h);  
  l <- list("A1A1h"=A1A1h,"A1A1H"=A1A1H,"A1A2h"=A1A2h,"A1A2H"=A1A2H,"A2A2h"=A2A2h,"A2A2H"=A2A2H) 			
  
  df <- data.frame( "gen"=(1:1440), "a370"=1/2/pop)
  
  for (i in 1:gen){ 
    l=ng(l$A1A1h,l$A1A1H,l$A1A2h,l$A1A2H,l$A2A2h,l$A2A2H,alpha=alpha,d=d,b0=b0,b1=b1,b2=b2,b3=b3,s=s,compen=compen,f=f)
    v370=l$A1A1h+l$A1A1H+(l$A1A2h+l$A1A2H)/2  #ancestral allele
    a370=1-v370
    df[(i+1440-gen),2]<-a370
  }
    print(a370)
return(df)}

ns0.074<-data.frame(f(b0=0,b1=0.5,b2=0.5,b3=1,d=0.5,s=0.074,compen=3.68,pop=80000,alpha=0,f=0,gen=1440))
ns0.066<-data.frame(f(b0=0,b1=0.5,b2=0.5,b3=1,d=0.5,s=0.066,compen=3.68,pop=80000,alpha=0,f=0,gen=1440))
ns0.066_f_alpha_1440<-data.frame(f(b0=0,b1=0.5,b2=0.5,b3=1,d=0.5,s=0.066,compen=3.68,pop=80000,alpha=0.005,f=0.056,gen=1440))
ns0_f_alpha_1440<-data.frame(f(b0=0,b1=0.5,b2=0.5,b3=1,d=0.5,s=0,compen=3.68,pop=80000,alpha=0.065,f=0.021,gen=1440))
ns0_f_alpha_1165<-data.frame(f(b0=0,b1=0.5,b2=0.5,b3=1,d=0.5,s=0,compen=3.68,pop=80000,alpha=0.0705,f=0.0342,gen=1165))
ns0.066_b_alpha_cb<-data.frame(f(b0=0.25,b1=0.40,b2=0.40,b3=1,d=0.5,s=0.066,compen=3.68,pop=80000,alpha=0.01,f=0,gen=1440))
ns0.066_b_alpha_cu<-data.frame(f(b0=0,b1=0,b2=0.945,b3=1,d=0.5,s=0.066,compen=3.68,pop=80000,alpha=0.005,f=0,gen=1440))
ns0_b_alpha<-data.frame(f(b0=0.14,b1=0.2,b2=0.27,b3=1,d=0.5,s=0,compen=3.68,pop=80000,alpha=0.052,f=0,gen=1440))

ns0.066_f_b_alpha<-data.frame(f(b0=0,b1=0.65,b2=0.65,b3=.95,d=0.5,s=0.066,compen=3.68,pop=80000,alpha=0.02,f=0.126,gen=1440))





plot<-data.frame(gen=-1439:0,ns0.074[,2],ns0.066[,2],
                 ns0_f_alpha_1440[,2],ns0_f_alpha_1165[,2],ns0_b_alpha[,2],
                ns0.066_f_alpha_1440[,2],ns0.066_b_alpha_cb[,2],ns0.066_b_alpha_cu[,2],
                ns0.066_f_b_alpha[,2])
colnames(plot) <- c("gen", LETTERS[1:9])

plot2 <- melt(plot,id='gen')  # convert to long format

cbPalette <- c("#000000", "#0072B2","#E69F00","darkolivegreen1", "#56B4E9", "#009E73", "#F0E442",  "#CC79A7", "#D55E00")

theme_set(theme_bw())

p<- ggplot(data=plot2,aes(x=gen, y=(value), colour=variable)) +   geom_line() +theme(legend.key=element_blank())
p<- p + scale_x_continuous(expand = c(0,0),breaks=c(-1439,-1000,-500,0)) +
         scale_y_continuous(breaks=c(0,0.250,0.500,0.750,0.916))
p<- p + geom_hline(y=0.916, colour = "red",linetype="dashed")
p <- p + ylab("EDAR370A  freqeucny")
p <- p+  scale_colour_manual(values=cbPalette,name = "Model")


plot3<- melt(plot[(1:240),],id='gen')
m1<- ggplot(data=plot3,aes(x=gen, y=value, colour=variable)) +   geom_line()
m1<-m1 + scale_x_continuous(expand = c(0,0),breaks=c(-1439,-1300,-1200))+ 
    scale_colour_manual(values=cbPalette) +
    theme(legend.position="none") +
    labs(x=NULL, y=NULL)

plot3<- melt(plot[(240:440),],id='gen')
m1.5<- ggplot(data=plot3,aes(x=gen, y=value, colour=variable)) +   geom_line()
m1.5<-m1.5 + scale_x_continuous(expand = c(0,0),breaks=c(-1200,-1100,-1000))+ 
  scale_colour_manual(values=cbPalette) +
  theme(legend.position="none") +
  labs(x=NULL, y=NULL)

plot3<- melt(plot[(440:640),],id='gen')
m2<- ggplot(data=plot3,aes(x=gen, y=value, colour=variable)) +   geom_line()
m2<-m2 + scale_x_continuous(expand = c(0,0),breaks=c(-1000,-900,-800))+ 
  scale_colour_manual(values=cbPalette) +
  theme(legend.position="none") +
  labs(x=NULL, y=NULL)

plot3<- melt(plot[(640:840),],id='gen')
m3<- ggplot(data=plot3,aes(x=gen, y=value, colour=variable)) +   geom_line()
m3<-m3 + scale_x_continuous(expand = c(0,0),breaks=c(-800,-700,-600))+ 
  scale_colour_manual(values=cbPalette) +
  theme(legend.position="none") +
  labs(x=NULL, y=NULL)

plot3<- melt(plot[(840:1040),],id='gen')
m4<- ggplot(data=plot3,aes(x=gen, y=value, colour=variable)) +   geom_line()
m4<-m4 +  scale_x_continuous(expand = c(0,0),breaks=c(-600,-500,-400))+
  scale_colour_manual(values=cbPalette) +
  theme(legend.position="none") +
  labs(x=NULL, y=NULL)

plot3<- melt(plot[(1040:1440),],id='gen')
m5<- ggplot(data=plot3,aes(x=gen, y=value, colour=variable)) +   geom_line()
m5<-m5 +  scale_x_continuous(expand = c(0,0),breaks=c(-400,-200,0))+
  scale_colour_manual(values=cbPalette) +
  theme(legend.position="none") +
  labs(x=NULL, y=NULL)

plot3<- melt(plot[(1240:1440),],id='gen')
m6<- ggplot(data=plot3,aes(x=gen, y=value, colour=variable)) +   geom_line()
m6<-m6 +  scale_x_continuous(expand = c(0,0),breaks=c(-200,-100,0))+
  scale_colour_manual(values=cbPalette) +
  theme(legend.position="none") +
  labs(x=NULL, y=NULL)



library(grid)


vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)


# One figure in row 1 and two figures in row 2
grid.newpage()
pushViewport(viewport(layout = grid.layout(5, 6)))
print(p, vp = vplayout(1:3, 1:6))
print(m1, vp = vplayout(4:5, 1))
print(m1.5, vp = vplayout(4:5, 2))
print(m2, vp = vplayout(4:5, 3))
print(m3, vp = vplayout(4:5, 4))
print(m4, vp = vplayout(4:5, 5))
print(m5, vp = vplayout(4:5, 6))
#print(m6, vp = vplayout(4:5, 6))

