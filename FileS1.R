# R Script to create a function that generates random individual-level data under a PT-SSH-like process and analyze them.

library(jtools)
library(truncnorm)
library(ggplot2)
library(ggeffects)
library(car)
library(interactions)
library(AICcmodavg)

# First of all, set your working directory. The function generates some figure that will be located in this working directory.
# setwd(D:/Dropbox/myfolder")


# Multiple plot function (this is just to paste different plots together)
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#########################################################################################################
# This is the function to generate and analyze random individual-level data under a PT-SSH-like process #
#########################################################################################################

# TRANSLATIONAL TABLE - The arguments refer to:
# - n.ses, the number of mating occasions, this may be more or less realistic depending on the system. It is to discretize the events and, therefore, to facilitate the simulation. 
# - max.entries, the maximum number of males that can enter into the breeding ground at a certain occasion.
# - min.indFit, the minimum number of grandoffspring a female can be predicted to obtain if she mated with the male with the lowest quality. In the manuscript, it corresponds to α in eq. 1 and eq. 3.
# - min.matingP, the minimum chance a male has of finding a mate at a certain mating occasion. In the manuscript, it corresponds to αp in eq. 2.
# - beta.indFit, the slope of the effect of male quality on females’ indirect benefits. In the manuscript, it corresponds to β in eq. 1 and eq. 3.
# - beta.mating, the slope of the effect of male quality on his chance of mating at a certain mating occasion. The greater its value the greater is the advantage in the probability of mating of a male with a given surplus in quality. It can be seen as a proxy of male-to-male competition. In the manuscript, it corresponds to βp in eq. 2.
# - C, this is the cost of mate sharing for the indirect benefits of a secondary female. This is also the cost a female is capable of predicting if she chooses an already-mated male. Choosing an already-mated male means that is 100% sure that she will be polygamously-mated. In the manuscript, it corresponds to C in eq. 1 and Csec in eq. 3. 
# - Cprim, this is the cost of mate sharing for the indirect benefits of a primary female. This cost acts on her indirect benefits but she cannot predict it at the moment of mating (as she mates with an unmated male). In the manuscript, it corresponds to Cprim in eq. 3.
# - sharingCosts, if =1 only secondary females pay a cost (-C) for mate sharing; if =2, the there is a cost (Cprim) also for primary females.  


PTSSH.fn<- function(n.ses=100,max.entries=100,min.indFit=2,min.matingP=0.2,beta.indFit=2,beta.mating=2,C=-0.4,sharingCosts=1,Cprim=-0.2){
  # Although this is a continuous process, to facilitate calculations the breeding sessions are discretized into "n.ses" choice sessions that is when a female chooses among the available males.
  if((sharingCosts==2)&(is.na(Cprim))) return("Error: if you are simulating sharingCosts=2 you need to set a value for Cprim")
  options(warn=-1)
  set.seed(223457)
  p.entries<- plogis(-1*(as.numeric(scale(1:n.ses)))-1*(as.numeric(scale(1:n.ses))^2))# This is the probability a male has of arriving to the breeding ground. To make it more realistic this probability has been set at its maximum value in the intermediate sessions.
  entries<- rbinom(n.ses,max.entries,p.entries)+1 # This is the session-specific numbers of individuals that enter into the population. These numbers are obtained by combining the individual probabilities of entering (p.entries) into the population and the number of individuals available to do it (max.entries).
  males<- matrix(0,sum(entries),n.ses) # Creating an empty matrix where the breeding histories of all the males will be stored.
  alpha.indFit<- log(min.indFit) # log-scale intercept for predicted female fitness
  alpha.mating<- qlogis(min.matingP) # logit-scale intercept for male probability of mating
  set.seed(223457)# this is useful for creating simulations or random objects that can be reproduced.
  entries<- rbinom(n.ses,max.entries,p.entries)+1 # This is the session-specific numbers of individuals that enter into the population. These numbers are obtained by combining the individual probabilities of entering (p.entries) into the population and the number of individuals available to do it (max.entries).
  males<- matrix(0,sum(entries),n.ses) # Creating an empty matrix where the breeding histories of all the males will be stored.
  qual<- rtruncnorm(n=sum(entries), a=(-C+0.01), b=(-C+1), mean=(-C+0.5), sd=0.1) # this is the vector with male qualities. Note that for the sake of simplicity I have not defined any correlation between the date at which an individual arrives to the breeding ground and its quality. However, this is just a baseline scenario, more complexity (and realism) may be added by including time-varying environmental covariates as predictors of the male mating probability and female indirect benefits.
  # Next steps are to populate the males matrix with "1"s representing their arrival to the breeding ground.
  entries0<- c(1,cumsum(entries)+1)
  entries1<- cumsum(entries)
  for(i in 1:(length(entries))){
    males[entries0[i]:entries1[i],i]<- 1
  }
  # In the next steps the males matrix is filled with the breeding histories of each individual. 
  for(j in 2:length(entries)){
    qualMon<-  qual[males[,j-1]==1] # at a certain session, a vector with the quality of all the unmated males available for female choice
    qualPol<-  qual[males[,j-1]==2] # at a certain session, a vector with the quality of all the already-mated males available for female choice.
    # Next steps are to create a single dataset (avMales) that contains all the males in the breeding ground that are available to mate.
    mon<- data.frame(qual=qualMon,ms=rep("mon",length(qualMon)))
    pol<- data.frame(qual=qualPol,ms=rep("pol",length(qualPol)))
    avMales<- rbind(mon,pol)
    msPenal<- ifelse(avMales$ms=="pol",C,0) # note that the predictable cost for female choice applies only to a male that is already-mated 
    indFit<- exp(alpha.indFit+beta.indFit*avMales$qual+msPenal)
    indFitOrd<- vector()
    indFitOrd<- (sort(indFit,decreasing=T)) # this represents a vector with the values of female fitness from the worse to the best available-to-mate male (according to their quality and mating status). 
    minIndFit<- numeric()
    minIndFit<-  indFitOrd[entries[j-1]] # This is the minimum quality a male needs to find a mate (at the net of its mating status - it will be clear in the next steps). Here, I am considering that the number of females is a limiting factor and, for simplicity, it is exactly the same of new males that arrive to the breeding ground. All the new females find a mate but only a fraction of new males find a mate. That fraction depends on the available to mate male in that session that is the sum of new males and those that did not find a mate in the previous sessions. This seems a realistic approximation as it entails that, over time, the male-to-male competition increases because it increases the number of competitors. In the next loop the male matrix is filled with: 2= already-mated male, and 3= polygynous male (only the case of max 2 females per male is contemplated).
    for (z in 1:nrow(males)){
      if(males[z,j-1]==0) {next} # if this male has just arrived, skip to the next male (next row in the matrix). This is because the initial state of each male when it enters the population must be a 1, i.e. as an unmated male.
      if(males[z,j-1]==1) { # this is the case when a male was an unmated male in the previous breeding session
        indFitz<- numeric()
        indFitz<- exp(alpha.indFit+beta.indFit*qual[z]) # therefore, the expected fitness for a female mating with this male is proportional to its quality (on the log scale) without any handicap related to the polygamy status (because it is still an unmated male)
        ifelse(indFitz>=minIndFit,males[z,j]<-rbinom(1,1,p=plogis(alpha.mating+beta.mating*indFitz))+1,males[z,j]<-1) # the male will have a chance to mate only if the expected fitness is superior to the minimum male quality specific for this breeding session. To calculate this value, I have considered that the number of females corresponds to the number of males that are new entries, say "n", and that only the first n males with the highest quality have a chance to mate.  
      }
      if(males[z,j-1]==2) { # this is the case of already-mated males. Their chance of finding a mate will follow the same competition but they have the handicap (from the female choice standpoint) of being already-mated.
        indFitz<- numeric()
        indFitz<- exp(alpha.indFit+beta.indFit*qual[z]+C)
        ifelse(indFitz>=minIndFit,males[z,j]<-rbinom(1,1,p=plogis(alpha.mating+beta.mating*indFitz))+2,males[z,j]<-2)
      }
      if(males[z,j-1]==3) { # this is the case of polygynous males. They have no chance of finding a further mate.
        males[z,j]<- males[z,j-1]
      }
    }
  }
  
  # Next steps are to create a dataset that will contain the breeding histories of females. As it happens in most cases when one deals with observational data (e.g. data from nest-boxes), these data inform uniquely on the breeding histories of males and females that have mated at some time during the breeding season (not of those who did not find a mate).
  arrivDate<- vector()
  monDate<- vector()
  polDate<- vector()
  for(i in 1:nrow(males)){
    arrivDate[i]<- min(which(males[i,] > 0))
    if(is.infinite(min(which(males[i,] > 1)))==TRUE){monDate[i]<- NA}
    else{monDate[i]<- min(which(males[i,] > 1))}
    if(is.infinite(min(which(males[i,] > 2)))==TRUE){polDate[i]<- NA}
    else{polDate[i]<- min(which(males[i,] > 2))}
  }
  # Translating the males histories into those of females:
  males<- cbind(males,qual,arrivDate,monDate,polDate)
  malesData<- as.data.frame(males)
  malesMat<- malesData[!is.na(malesData$monDate),]
  malesMat$malID<- rownames(malesMat) # this will be useful to define male identity as a random intercept (polygynous males appear twice in the data set). 
  malesMat$ms<- ifelse(is.na(malesMat$polDate),"mon","pol")
  femMon<- malesMat[is.na(malesMat$polDate),]
  femMon$date<- femMon$monDate
  if(nrow(malesMat[!is.na(malesMat$polDate),])==0) return("Error: Maybe the C is too big or the beta.indFit too small to make the polygyny threshold possible. If the effect size of quality (beta.indFit) on the female indirect fitness is small, you need a small polygyny cost (C) value to make the female choice of an already-mated male possible to happen.")
  femPol<- malesMat[!is.na(malesMat$polDate),]
  femPrim<- femPol
  femSec<- femPol
  femPrim$date<- femPrim$monDate
  femPrim$ms<- "prim"
  femSec$date<- femSec$polDate
  femSec$ms<- "sec"
  fem<- rbind(femMon,femPrim,femSec)
  fem$ms<-as.factor(fem$ms)
  fem<- fem[,(ncol(fem)-7):ncol(fem)]
  fem<- fem[,c(1,2,6,7,8)]
  tablems<- table(fem$ms)
  malesMat$ms<- as.factor(malesMat$ms)
  
  pol<- malesMat[malesMat$ms=="pol",]
  pol<- pol[order(pol$qual),]
  pol$inds<- seq(0,1,length.out=nrow(pol))
  mon<- malesMat[malesMat$ms=="mon",]
  mon<- mon[order(mon$qual),]
  mon$inds<- seq(0,1,length.out=nrow(mon))
  minQP<- min(pol$qual)
  nPol<- nrow(pol)
  
  fem$monpol<- as.factor(ifelse(fem$ms=="sec","pol","mon"))
  # str(fem)
  femMon<- fem[fem$monpol=="mon",]
  femMon<- femMon[order(femMon$qual),]
  femMon$inds<- seq(0,1,length.out=nrow(femMon))
  femPol<- fem[fem$monpol=="pol",]
  femPol<- femPol[order(femPol$qual),]
  femPol$inds<- seq(0,1,length.out=nrow(femPol))
  femMon1<- fem[fem$ms=="mon",]
  femMon1<- femMon1[order(femMon1$qual),]
  femMon1$inds<- seq(0,1,length.out=nrow(femMon1))
  femSec<- fem[fem$ms=="sec",]
  femSec<- femSec[order(femSec$qual),]
  femSec$inds<- seq(0,1,length.out=nrow(femSec))
  
  ####################################
  #                                  #
  #  Fig. S1 - Male quality overlap  #
  #                                  #
  ####################################
  
  tiff("FigS1.tiff", width = 16.6, height = 16.6, units = 'cm', res = 600, compression="lzw")
  # png("Fig.S1.png",width=600,height=960)
  par(mfrow=c(3,1),mar=c(4,8,2,2))
  plot(mon$inds,mon$qual,pch=19,col=adjustcolor("black",0.5),main="",ylim=range(malesMat$qual),ylab="Male  Quality",xlab="Males",xaxt="n",bty="l",cex.lab=1)
  points(pol$inds,pol$qual,pch=19,col=adjustcolor("red",0.5),main="")
  abline(h=minQP,col=adjustcolor("red",0.5),lwd=.8)
  legend("topleft",legend=c("Polygynous male","Monogamous male"),col=c(adjustcolor("red",0.5),adjustcolor("black",0.5)),pch=c(19,19),cex=1)
  
  plot(femMon$inds,femMon$qual,pch=19,col=adjustcolor("black",0.5),main="",ylim=range(fem$qual),ylab="Male  Quality",xlab="Females",xaxt="n",bty="l",cex.lab=1)
  points(femPol$inds,femPol$qual,pch=19,col=adjustcolor("red",0.5),main="")
  abline(h=min(femPol$qual),col=adjustcolor("red",0.5),lwd=.8)
  legend("topleft",legend=c("Mating with already-mated male","Mating with unmated male"),col=c(adjustcolor("red",0.5),adjustcolor("black",0.5)),pch=c(19,19),cex=1)
  
  plot(femMon1$inds,femMon1$qual,pch=19,col=adjustcolor("black",0.5),main="",ylim=range(fem$qual),ylab="Male  Quality",xlab="Females",xaxt="n",bty="l",cex.lab=1)
  points(femSec$inds,femSec$qual,pch=19,col=adjustcolor("red",0.5),main="")
  abline(h=min(femSec$qual),col=adjustcolor("red",0.5),lwd=.8)
  legend("topleft",legend=c("Polygamously-mated female (prim. and sec.)","Monogamous female"),col=c(adjustcolor("red",0.5),adjustcolor("black",0.5)),pch=c(19,19),cex=1)
  dev.off()
  
  
  if(sharingCosts==1){
    #########################################################################
    # sharingCosts  1 - only secondary females pay a cost in direct fitness #
    #########################################################################
    # the number of grandoffspring of a female depends on male quality and male mating status in such a way that only secondary females are penalized by the mate polygynous status. For instance, with species with biparental care this scenario holds when a polygynous male takes care only of the offspring from the primary female.
    fem1<- fem
    fem1$date<- as.numeric(scale(fem1$date))
    fem1$msPenal<- ifelse(fem1$ms=="sec",C,0)
    fem1$indFit<- rpois(nrow(fem1),exp(alpha.indFit+beta.indFit*fem1$qual+fem1$msPenal))
    fem1$monpol<- as.factor(ifelse(fem1$ms=="sec","pol","mon"))
    # Next, we build two models to explain the variation in fem1ale Indirect  fitness by considering as predictors: 1) the male mating status and, 2) male quality and mating status. The 1) has often been done in observational studies to test the sexy son hypothesis. The expectation is that secondary (and sometimes the primary, see the ms text) fem1ales have an equal or superior Indirect  fitness than monogamous fem1ales. The 2) is what I recommend to do to evaluate the PT-SSH, however this model serves to evaluate compatible expectations generated by the PT-SSH. As explained in the text, a critical test can be achieved only by an experimental approach. 
    GLMqual<- glm(indFit~qual,data=fem1,"poisson")
    GLMqualms<- glm(indFit~qual+ms,data=fem1,"poisson")
    GLMms<- glm(indFit~ms,data=fem1,"poisson")
    GLMqualmonpol<- glm(indFit~qual+monpol,data=fem1,"poisson")
    GLMmonpol<- glm(indFit~monpol,data=fem1,"poisson")
    
    Cand.models<- list()
    Cand.models[[1]]<- GLMqual
    Cand.models[[2]]<- GLMqualms
    Cand.models[[3]]<- GLMms
    Cand.models[[4]]<- GLMqualmonpol
    Cand.models[[5]]<- GLMmonpol
    
    Modnames<- paste(c("fitness ~ qual","fitness ~ qual + ms(mon,prim,sec)","fitness ~ ms(mon,prim,sec)","fitness ~ qual + ms(mon,pol)","fitness ~ ms(mon,pol)"),sep=" ")
    
    print(AIC.TabOrd<- aictab(cand.set=Cand.models,modnames=Modnames))
    
    PT<- deltaMethod(GLMqualmonpol, "b2/b1", parameterNames= paste("b", 0:2, sep=""))# this is the estimate of the polygyny threshold based on the ratio between the fitness cost for the female choosing an already-mated male (with respect to one choosing an unmated male of identical quality) and the effect size of quality on female fitness. The uncertainty on this estimate is calculated by using the delta method.
    
    # For Fig. S2:
    trueMon<- exp(log(min.indFit)+beta.indFit*mean(fem1$qual))
    truePrim<- exp(log(min.indFit)+beta.indFit*mean(fem1$qual))
    trueSec<- exp(log(min.indFit)+beta.indFit*mean(fem1$qual)+C)
    datims<- data.frame(ms=1:3,indFitT=c(trueMon,truePrim,trueSec))
    datimonpol<- data.frame(monpol=1:2,indFitT=c(trueMon,trueSec))
    
    output<<- list(GLMqualms=GLMqualms,GLMms=GLMms,GLMqualmonpol=GLMqualmonpol,GLMmonpol=GLMmonpol,AIC.TabOrd=AIC.TabOrd,males=males,tablems=tablems,PT=PT,malesMat=malesMat,data=data,fem=fem)}
  
  
  if(sharingCosts==2){
    #########################################################################
    # sharingCosts  2 - prim. and sec. females pay a cost in direct fitness #
    #########################################################################
    # the number of grandoffspring of a female depends on male quality and male mating status in such a way that secondary and primary females are penalized by the mate polygynous status. As it is often the case in bird species, secondary females pay a higher cost than primary females.
    fem2<- fem
    fem2$date<- as.numeric(scale(fem2$date))
    fem2$msPenal<- ifelse(fem2$ms=="sec",C,ifelse(fem2$ms=="prim",Cprim,0))
    fem2$indFit<- rpois(nrow(fem2),exp(alpha.indFit+beta.indFit*fem2$qual+fem2$msPenal))
    fem2$monpol<- as.factor(ifelse(fem2$ms=="sec","pol","mon"))
    fem2$monpol2<- as.factor(ifelse(fem2$ms=="mon","mon","pol"))
    # Next, we build two models to explain the variation in female indirect benefits by considering as predictors: 1) the male mating status and, 2) male quality and mating status. The 1) has often been done in observational studies to test the sexy son hypothesis. The expectation is that secondary (and sometimes the primary, see the ms text) fem2ales have an equal or superior indirect benefits than monogamous fem2ales. The 2) is what I recommend to do to evaluate the PT-SSH, however this model serves to evaluate compatible expectations generated by the PT-SSH. As explained in the text, a critical test can be achieved only by an experimental approach. 
    GLMqual<- glm(indFit~qual,data=fem2,"poisson")
    GLMqualms<- glm(indFit~qual+ms,data=fem2,"poisson")
    GLMms<- glm(indFit~ms,data=fem2,"poisson")
    GLMqualmonpol<- glm(indFit~qual+monpol,data=fem2,"poisson")
    GLMmonpol<- glm(indFit~monpol,data=fem2,"poisson")
    
    Cand.models<- list()
    Cand.models[[1]]<- GLMqual
    Cand.models[[2]]<- GLMqualms
    Cand.models[[3]]<- GLMms
    Cand.models[[4]]<- GLMqualmonpol
    Cand.models[[5]]<- GLMmonpol
    
    Modnames<- paste(c("fitness ~ qual","fitness ~ qual + ms(mon,prim,sec)","fitness ~ ms(mon,prim,sec)","fitness ~ qual + ms(mon,pol)","fitness ~ ms(mon,pol)"),sep=" ")
    
    print(AIC.TabOrd<- aictab(cand.set=Cand.models,modnames=Modnames))
    
    PT<- deltaMethod(GLMqualmonpol, "b2/b1", parameterNames= paste("b", 0:2, sep=""))# this is the estimate of the polygyny threshold based on the ratio between the fitness cost for the female choosing an already-mated male (with respect to one choosing an unmated male of identical quality) and the effect size of quality on female fitness. The uncertainty on this estimate is calculated by using the delta method.
    
    # For Fig. S2:
    trueMon<- exp(log(min.indFit)+beta.indFit*mean(fem2$qual))
    truePrim<- exp(log(min.indFit)+beta.indFit*mean(fem2$qual)+Cprim)
    trueSec<- exp(log(min.indFit)+beta.indFit*mean(fem2$qual)+C)
    datims<- data.frame(ms=1:3,indFitT=c(trueMon,truePrim,trueSec))
    datimonpol<- data.frame(monpol=1:2,indFitT=c(trueMon,trueSec))
    
    output<<- list(GLMqualms=GLMqualms,GLMms=GLMms,GLMqualmonpol=GLMqualmonpol,GLMmonpol=GLMmonpol,AIC.TabOrd=AIC.TabOrd,males=males,tablems=tablems,PT=PT,malesMat=malesMat,data=data,fem=fem)}
  
  ####################################
  #                                  #
  #   Fig. S2 - Models' predictions  #
  #                                  #
  ####################################
  
  p1<- plot(ggpredict(GLMms,"ms"))+geom_point(size=2)+
    geom_point(data= datims, aes(x = ms+0.01, y = indFitT), colour= 'red', size = 2)+
    theme_gray()+
    ggtitle("Model {fitness ~ ms(mon,prim,sec)}")+
    labs(x="",y=expression("Ind.  benefits  "["(no.  grandoff.)"]))+
    theme(plot.title = element_text(face="bold", size=12, hjust=0))+
    theme(axis.title = element_text(size=10))+
    theme(axis.text = element_text(face="bold", size=10)) 
  
  p2<- plot(ggpredict(GLMqualms,"ms"))+geom_point(size=2)+
    geom_point(data= datims, aes(x = ms+0.01, y = indFitT), colour= 'red', size = 2)+
    theme_gray()+
    ggtitle("Model {fitness ~ ms(mon,prim,sec) + qual}")+
    labs(x="",y=expression("Ind.  benefits  "["(no.  grandoff.)"]))+
    theme(plot.title = element_text(face="bold", size=12, hjust=0))+
    theme(axis.title = element_text(size=10))+
    theme(axis.text = element_text(face="bold", size=10)) 
  
  p3<- plot(ggpredict(GLMmonpol,"monpol"))+geom_point(size=2)+
    geom_point(data= datimonpol, aes(x = monpol+0.01, y = indFitT), colour= 'red', size = 2)+
    theme_gray()+
    ggtitle("Model {fitness ~ ms(mon,pol)}")+
    labs(x="",y=expression("Ind.  benefits  "["(no.  grandoff.)"]))+
    theme(plot.title = element_text(face="bold", size=12, hjust=0))+
    theme(axis.title = element_text(size=10))+
    theme(axis.text = element_text(face="bold", size=10)) 
  
  p4<- plot(ggpredict(GLMqualmonpol,c("monpol")))+geom_point(size=2)+
    geom_point(data= datimonpol, aes(x = monpol+0.01, y = indFitT), colour= 'red', size = 2)+
    ggtitle("Model {fitness ~ ms(mon,pol) + qual}")+
    theme_gray()+
    labs(x="",y=expression("Ind.  benefits  "["(no.  grandoff.)"]))+
    theme(plot.title = element_text(face="bold", size=12, hjust=0))+
    theme(axis.title = element_text(size=10))+
    theme(axis.text = element_text(face="bold", size=10)) 
  
  p5<- interact_plot(GLMqualms,modx=ms,pred=qual,legend.main="Female mating status",x.label="Male Quality",y.label=expression("Ind.  benefits  "["(link  scale)"]),modx.labels=c("Monogamous","Primary","Secondary"),outcome.scale="link")
  p5<- p5+theme(axis.title.y = element_text(size = rel(1), angle = 90))
  p5<- p5 + theme(axis.title.x = element_text(size = rel(1), angle = 00))
  p5<- p5 + theme(legend.text=element_text(size=rel(1)))
  
  p6<- interact_plot(GLMqualmonpol,modx=monpol,pred=qual,legend.main="Female mating with:",x.label="Male Quality",y.label=expression("Ind.  benefits  "["(link  scale)"]),modx.labels=c("Unmated","Already-mated"),outcome.scale="link")
  p6<- p6+theme(axis.title.y = element_text(size = rel(1), angle = 90))
  p6<- p6 + theme(axis.title.x = element_text(size = rel(1), angle = 00))
  p6<- p6 + theme(legend.text=element_text(size=rel(1)))
  newdata<- GLMqualmonpol$data
  span.x<- seq(min(newdata$qual),max(newdata$qual),0.001)
  span.y<- seq(min(log(GLMqualmonpol$fitted.values)),max(log(GLMqualmonpol$fitted.values)),0.001)
  newdataPol<- newdata[newdata$monpol=="pol",]
  x2<- newdataPol[which.min(abs(newdataPol$qual-quantile(span.x,0.5))),]$qual
  x1<- x2-abs(PT$Estimate)
  y1<- y2<- predict(GLMqualmonpol,newdata=newdata[newdata$qual==x2&newdata$monpol=="pol",])
  df <- data.frame(x1=x1, x2=x2, y1=y1, y2=y2)
  p6<- p6 + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),inherit.aes=FALSE,data = df)
  p6<- p6 + annotate("text", x= quantile(span.x,0.75), y= y1, label= paste(("PT = "), round(abs(PT$Estimate),2),sep=" "),size=3)
  
  tiff("FigS2.tiff", width = 26, height = 15, units = 'cm', res = 600, compression="lzw")
  # png("FigS2.png",width=1280,height=800)
  multiplot(p1,p3,p5,p2,p4,p6,cols=2)
  dev.off()
  
  data<- GLMmonpol$data
  
  ####################################
  #                                  #
  #    Fig. S3 - Females' benefits   #
  #                                  #
  ####################################
  
  tiff("FigS3.tiff", width = 16.6, height = 16.6, units = 'cm', res = 600, compression="lzw")
  # png("FigS3.png",width=1280*0.8,height=800*0.8)
  par(mar=c(4,8,2,2))
  boxplot(data$indFit~data$ms,xlab="Female mating status",ylab=expression("Ind.  benefits  "["(no.  grandoff.)"]),cex.lab=1.5,cex.axis=1.5)
  abline(h=mean(data[data$ms=="mon",]$indFit),lty=2)
  dev.off()
  
  invisible(output)
  options(warn=0)
}


################################################
#     Now you can use the PTSSH.fn function    #
################################################
# Use it to generate random individual data from a PT-SSH process and running models to investigate how the female indirect fitness depends on the interplay between her mating status and the quality of her mate.  

PTSSH.fn(n.ses=50,max.entries=200,min.indFit=2,min.matingP=0.5,beta.indFit=2,beta.mating=2,C=-0.5,Cprim=-0.1,sharingCosts=2)
# The first result you can look at is an AIC table ranking five different models of female indirect fitness dependin on: (i) male quality (GLMqual); (ii) female mating status (mon,prim,sec) (GLMms); (iii) female choice of an unmated or an already-mated male (mon+prim,sec) (GLMmonpol); (iv) male quality + female mating status (mon,prim,sec) (GLMqualms); (v) male quality + female choice of an unmated or an already-mated male (mon+prim,sec) (GLMqualmonpol). The lower the AIC, the more support the model yields.

list2env(output, .GlobalEnv)# this is to volcate all the objects created by PTSSH.fn in your workspace

PT # This is the estimate of polygyny treshold, i.e. the difference in quality an already-mated male must have with respect to a competitor unmated male to have a higher chance of being chosen by a female. 
tablems # This is a table with the number of monogamous (mon), primary (prim) and secondary (sec) females present in the simulated dataset.
