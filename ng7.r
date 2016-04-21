model<-function(pop=80000,yr=36000,p1=0.916,d=0.5,s=0.066,sk=0.074,pref.from.beginning= TRUE){

	#pop=starting population, yr=time since mutation, p1 = known current mutation frequency
	#d=dominance parameter of mutation (semi-dominant=0.5)
	#s=fitting selection coefficient (pure nature selection frequency,assumed);
	#sk= known selection coefficient (may include cultural sexual pressure, thus higher than pure nature selection coefficient)
	#pref.from.beginning = TRUE assumes the emergence of H preference coincide the emergence of A2, else pref emerge from a later generation 
	pt=1/2/pop#starting allele frequency for mutation
	gen=yr/25#assuming generation time =25yr
	h=1-1/pop #initial frequency of people who does not have a preference of mutation phenotype
	#evo1 calculates change in allele frequency for next generation, assuming only natural selection
	evo1 <- function (pt, sk=sk, d=d, compen) { 
			qt=1-pt
			sx=sk/compen #adjusted selection coefficient used for calculation
			p <- (pt^2*(1+sx)+pt*qt*(1+d*sx))/(1+sx*(pt^2+2*d*pt*qt))
			return(p)
		}
	#evo2 calculates "compensation factor" for selection pressure 
	evo2<- function(pop,yr,p1){ 
		ptx=0 #ptx= the p(allele frequency) produced by the model closest to p1 #initialisingng
		for (compen in seq(1,5, by=0.01)){ #fitting values for compensation factor in range 1-5
			pt=1/2/pop
			for (i in 1:gen){
				pt<- evo1(pt=pt,sk=sk,d=d,compen=compen)}
			if (abs(pt-p1)<abs(ptx-p1)){
				cmpx<-compen
				ptx<-pt
			}
		}
		print("compensation factor")
		return(cmpx)
	}
	#ng calculates pheno-genotype in next generation
	ng <- function(x1,x2,x3,x4,x5,x6,alpha,d,b0,b1,b2,b3,s,compen,f){ #s=true selection coefficient comp=compensation factor
		sx<-s/compen
		z<-1+alpha*d*x3+alpha*d*x4+alpha*x5+alpha*x6; #normalising factor
		a<-(1+alpha)
		da<-(1+d*alpha)
		#next gen pheno-genotype, without selection
		A1A1h<-b3*(x1*x1+x1*x3+1/4*x3*x3)+b2*(x1*x2+0.5*x2*x3+0.5*x1*x4+0.25*x3*x4)+b1/z*(x2*x1+0.5*x2*x3*da+0.5*x1*x4+0.25*da*x3*x4)+b0/z*(x2*x2+0.5*x2*x4*da+0.5*x4*x2+0.25*x4*x4*da)
		A1A1H<-(1-b3)*(x1*x1+x1*x3+1/4*x3*x3)+(1-b2)*(x1*x2+0.5*x2*x3+0.5*x1*x4+0.25*x3*x4)+(1-b1)/z*(x2*x1+0.5*x2*x3*da+0.5*x1*x4+0.25*da*x3*x4)+(1-b0)/z*(x2*x2+0.5*x2*x4*da+0.5*x4*x2+0.25*x4*x4*da)
		A1A2h<-b3*(2*x1*x5+x1*x3+x3*x5+0.5*x3*x3)+b2*(x1*x6+x2*x5+0.5*(x1*x4+x2*x3+x3*x4+x3*x6+x4*x5))+b1/z*(x2*x5*a+x1*x6+0.5*(x2*x3*da+x1*x4+x3*x4*da+x4*x5*a+x3*x6*da))+b0/z*(x2*x6*a+x2*x6+0.5*(x2*x4*da+x2*x4+x4*x4*da+x4*x6*a+x6*x4*da))
		A1A2H<-(1-b3)*(2*x1*x5+x1*x3+x3*x5+.5*x3*x3)+(1-b2)*(x1*x6+x2*x5+0.5*(x1*x4+x2*x3+x3*x4+x3*x6+x4*x5))+(1-b1)/z*(x2*x5*a+x1*x6+0.5*(x2*x3*da+x1*x4+x3*x4*da+x4*x5*a+x3*x6*da))+(1-b0)/z*(x2*x6*a+x2*x6+0.5*(x2*x4*da+x2*x4+x4*x4*da+x4*x6*a+x6*x4*da))
		A2A2h<-b3*(x5*x5+x5*x3+1/4*x3*x3)+b2*(x5*x6+0.5*x6*x3+0.5*x5*x4+0.25*x3*x4)+b1/z*(x5*x6*a+0.5*x4*x5*a+0.5*x3*x6*da+0.25*da*x3*x4)+b0/z*(x6*x6*a+0.5*x6*x4*a+0.5*x6*x4*da+0.25*x4*x4*da)
		A2A2H<-(1-b3)*(x5*x5+x5*x3+1/4*x3*x3)+(1-b2)*(x5*x6+0.5*x6*x3+0.5*x5*x4+0.25*x3*x4)+(1-b1)/z*(x5*x6*a+0.5*x4*x5*a+0.5*x3*x6*da+0.25*da*x3*x4)+(1-b0)/z*(x6*x6*a+0.5*x6*x4*a+0.5*x6*x4*da+0.25*x4*x4*da)  
		#adding selection coefficient
		nor<-((1+sx)*(A2A2h+A2A2H)+(1+d*sx)*(A1A2h+A1A2H)+(A1A1h+A1A1H))#normalising
		A1A1hx<-A1A1h/nor
		A1A1Hx<-A1A1H/nor
		A1A2hx<-A1A2h*(1+d*sx)/nor
		A1A2Hx<-A1A2H*(1+d*sx)/nor
		A2A2hx<-A2A2h*(1+sx)/nor
		A2A2Hx<-A2A2H*(1+sx)/nor
		
		#adding oblique transmission
		A1A1H <- A1A1Hx+f*A1A1hx*(A1A1Hx+A1A2Hx+A2A2Hx)
		A1A2H <- A1A2Hx+f*A1A2hx*(A1A1Hx+A1A2Hx+A2A2Hx)
		A2A2H <- A2A2Hx+f*A2A2hx*(A1A1Hx+A1A2Hx+A2A2Hx)
		A1A1h <- A1A1hx*(1-f*(A1A1Hx+A1A2Hx+A2A2Hx))
		A1A2h <- A1A2hx*(1-f*(A1A1Hx+A1A2Hx+A2A2Hx))
		A2A2h <- A2A2hx*(1-f*(A1A1Hx+A1A2Hx+A2A2Hx))
		
		return(list("A1A1h"=A1A1h,"A1A1H"=A1A1H,"A1A2h"=A1A2h,"A1A2H"=A1A2H,"A2A2h"=A2A2h,"A2A2H"=A2A2H))
	}
	
	#apply evo2 to calculate compensation factor
	compen<-evo2(pop,yr,p1)
	print(compen)
	#starting genotype frequencies	
	A2=pt  #starting frequencies of mutation allele 
	A1=1-A2
	A1A1=A1*A1
	A1A2=A1*A2*2
	A2A2=A2*A2
	#allele frequencies today
	a2=p1
	a1=1-a2
	#fitting parameter range
	##fitting
	diff=1 #initialising least difference
	A1A1h=A1A1*h;
	A1A1H=A1A1*(1-h);
	A1A2h=A1A2*h;
	A1A2H=A1A2*(1-h);
	A2A2h=A2A2*h;
	A2A2H=A2A2*(1-h);
	for (b0 in #range){
		for (b1 in #range)){
			for (b2 in #range)){
				for (b3 in #range)){
					#if (b0>b1&b1<=b2){
					#	b0=b1}
					#if (b0>b2&b1>b2){
					#	b0=b2}
					#if (b3<b1){
					#	b3=b1}
					#if (b3<b2){
					#	b3=b2}	
						for (alpha in #range){
						for (f in #range){
							for (gen in seq(as.integer(yr/25),as.integer(0.8*(yr/25)),by=-25)){
								if (pref.from.beginning){
									gen <- as.integer(yr/25)
								}
								l <- list("A1A1h"=A1A1h,"A1A1H"=A1A1H,"A1A2h"=A1A2h,"A1A2H"=A1A2H,"A2A2h"=A2A2h,"A2A2H"=A2A2H)			
								for (i in 1:gen){ 
									l=ng(l$A1A1h,l$A1A1H,l$A1A2h,l$A1A2H,l$A2A2h,l$A2A2H,alpha=alpha,d=d,b0=b0,b1=b1,b2=b2,b3=b3,s=s,compen=compen,f=f)                
								}
								V370=l$A1A1h+l$A1A1H+(l$A1A2h+l$A1A2H)/2  #ancestral allele
								if (abs(V370- a1) < diff){
									diff <- abs(V370- a1)
									genx <- gen									
									pa<- c("f"=f,"b0"=b0,"b1"=b1,"b2"=b2,"b3"=b3,"gen"=genx,"alpha"=alpha)
									prefh=l$A1A1h+l$A1A2h+l$A2A2h
									prefH=1-prefh
									v370=l$A1A1h+l$A1A1H+(l$A1A2h+l$A1A2H)/2 
									a370=1-v370
									print (pa)
									print(c("h"=prefh,"f"=f,"ancestral allele"=v370,"gen"=genx, "diff"=diff))
								}
							}  
						}
					}
				}
			}
		}			
	}
	print("best parameter:")
	print(pa)
	print((c("h"=prefh,"H"=prefH,"ancestral allele"=v370,"derived allele"=a370)))
}

model(pop=80000,yr=36000,p1=0.916,d=0.5,s=0.066,sk=0.074,pref.from.beginning = TRUE)