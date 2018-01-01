pdf(file='GROWTH_10005_100_10_0.1_0.01.pdf',width=4.5,height=4.5);
gstable=read.table('GROWTH_10005_100_10_0.1_0.01.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("G0913","G0753","G0562","G0604","G0792","G0308","G0620","G0884","G0113","G0295")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='H1,H2,H3,H4_vs_L1,L2,L3,L4 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(16.3757586359,15.5082430436,15.3685458228,8.49397981418,7.54120649218,6.38007579132,8.19480548684,5.36013172842),c(8.18787931794,9.4398001135,6.32822475055,8.49397981418,6.03296519374,5.58256631741,2.04870137171,2.68006586421),c(9.0976436866,5.39417149343,2.71209632167,8.49397981418,4.52472389531,3.19003789566,1.02435068586,0.0),c(4.5488218433,7.41698580346,7.23225685778,4.24698990709,6.03296519374,3.98754736958,3.07305205757,0.0),c(13.6464655299,12.1368858602,10.8483852867,10.6174747677,0.0,2.39252842175,2.04870137171,2.68006586421),c(6.36835058062,7.41698580346,6.32822475055,2.12349495355,4.52472389531,4.78505684349,6.14610411513,0.0),c(0.0,2.02281431004,2.71209632167,2.12349495355,0.0,0.797509473916,2.04870137171,2.68006586421),c(13.6464655299,14.1597001702,16.27257793,23.358444489,0.0,7.17758526524,8.19480548684,8.04019759263),c(16.3757586359,8.09125724014,16.27257793,6.37048486064,15.0824129844,12.7601515826,13.3165589161,10.7202634568),c(2.72929310598,2.02281431004,4.52016053611,4.24698990709,0.0,0.797509473916,2.04870137171,0.0))
targetgene="G0913"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3.63905747464,7.41698580346,9.04032107222,10.6174747677,3.01648259687,2.39252842175,3.07305205757,2.68006586421),c(4.5488218433,2.69708574671,5.42419264333,4.24698990709,4.52472389531,0.0,2.04870137171,2.68006586421),c(0.90976436866,0.674271436679,0.904032107222,0.0,1.50824129844,0.0,1.02435068586,0.0),c(0.90976436866,1.34854287336,4.52016053611,4.24698990709,3.01648259687,0.0,3.07305205757,0.0),c(10.9171724239,9.4398001135,7.23225685778,8.49397981418,1.50824129844,4.78505684349,6.14610411513,0.0),c(6.36835058062,8.76552867682,3.61612842889,4.24698990709,0.0,4.78505684349,0.0,0.0),c(1.81952873732,2.02281431004,3.61612842889,8.49397981418,1.50824129844,0.0,1.02435068586,2.68006586421),c(2.72929310598,4.71990005675,5.42419264333,4.24698990709,1.50824129844,0.797509473916,1.02435068586,2.68006586421),c(3.63905747464,4.04562862007,4.52016053611,2.12349495355,0.0,1.59501894783,1.02435068586,2.68006586421),c(0.90976436866,2.02281431004,3.61612842889,2.12349495355,0.0,0.797509473916,1.02435068586,2.68006586421))
targetgene="G0753"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(20.9245804792,20.2281431004,22.6008026805,10.6174747677,18.0988955812,11.1651326348,13.3165589161,8.04019759263),c(1.81952873732,0.674271436679,6.32822475055,4.24698990709,0.0,0.0,0.0,5.36013172842),c(0.90976436866,1.34854287336,0.904032107222,2.12349495355,0.0,0.0,0.0,8.04019759263),c(1.81952873732,8.09125724014,9.04032107222,0.0,1.50824129844,3.98754736958,1.02435068586,5.36013172842),c(2.72929310598,2.69708574671,0.904032107222,0.0,0.0,0.797509473916,0.0,0.0),c(1.81952873732,3.37135718339,2.71209632167,2.12349495355,0.0,0.0,0.0,0.0),c(6.36835058062,10.7883429869,3.61612842889,2.12349495355,3.01648259687,2.39252842175,1.02435068586,8.04019759263),c(19.1050517419,14.8339716069,9.04032107222,10.6174747677,7.54120649218,3.98754736958,10.2435068586,8.04019759263),c(7.27811494928,8.76552867682,4.52016053611,10.6174747677,7.54120649218,1.59501894783,0.0,2.68006586421),c(3.63905747464,2.02281431004,0.904032107222,0.0,0.0,0.0,0.0,0.0))
targetgene="G0562"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(10.9171724239,10.7883429869,10.8483852867,12.7409697213,3.01648259687,3.19003789566,4.09740274342,0.0),c(6.36835058062,5.39417149343,1.80806421444,6.37048486064,0.0,2.39252842175,0.0,0.0),c(2.72929310598,3.37135718339,3.61612842889,6.37048486064,0.0,0.0,3.07305205757,2.68006586421),c(0.0,4.71990005675,5.42419264333,4.24698990709,1.50824129844,2.39252842175,0.0,0.0),c(0.0,0.0,2.71209632167,0.0,0.0,1.59501894783,0.0,0.0),c(1.81952873732,0.674271436679,0.904032107222,2.12349495355,0.0,1.59501894783,0.0,0.0),c(1.81952873732,2.02281431004,1.80806421444,6.37048486064,1.50824129844,0.0,1.02435068586,0.0),c(2.72929310598,0.674271436679,5.42419264333,2.12349495355,1.50824129844,1.59501894783,1.02435068586,0.0),c(0.90976436866,2.69708574671,0.0,0.0,0.0,0.797509473916,0.0,0.0),c(2.72929310598,0.0,2.71209632167,0.0,0.0,0.0,0.0,0.0))
targetgene="G0604"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(0.90976436866,0.674271436679,2.71209632167,2.12349495355,0.0,0.797509473916,2.04870137171,0.0),c(5.45858621196,4.04562862007,3.61612842889,8.49397981418,0.0,3.19003789566,2.04870137171,0.0),c(7.27811494928,6.06844293011,9.04032107222,2.12349495355,0.0,1.59501894783,4.09740274342,5.36013172842),c(19.1050517419,6.74271436679,2.71209632167,12.7409697213,12.0659303875,8.77260421307,8.19480548684,5.36013172842),c(10.0074080553,7.41698580346,7.23225685778,8.49397981418,1.50824129844,9.57011368699,5.12175342928,2.68006586421),c(8.18787931794,0.0,3.61612842889,2.12349495355,1.50824129844,0.0,1.02435068586,2.68006586421),c(4.5488218433,0.674271436679,1.80806421444,2.12349495355,0.0,0.797509473916,0.0,0.0),c(0.0,1.34854287336,0.0,0.0,0.0,0.0,0.0,0.0),c(4.5488218433,4.71990005675,8.136288965,8.49397981418,1.50824129844,0.0,3.07305205757,2.68006586421),c(3.63905747464,4.71990005675,3.61612842889,0.0,0.0,0.0,0.0,0.0))
targetgene="G0792"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(0.90976436866,6.74271436679,6.32822475055,2.12349495355,1.50824129844,0.797509473916,1.02435068586,2.68006586421),c(3.63905747464,6.06844293011,6.32822475055,4.24698990709,1.50824129844,2.39252842175,1.02435068586,0.0),c(0.90976436866,1.34854287336,0.0,6.37048486064,3.01648259687,1.59501894783,3.07305205757,5.36013172842),c(1.81952873732,8.09125724014,7.23225685778,8.49397981418,3.01648259687,8.77260421307,4.09740274342,0.0),c(4.5488218433,5.39417149343,1.80806421444,0.0,4.52472389531,0.797509473916,1.02435068586,5.36013172842),c(2.72929310598,0.674271436679,0.0,4.24698990709,3.01648259687,0.0,0.0,2.68006586421),c(0.0,0.674271436679,0.0,4.24698990709,0.0,0.0,0.0,0.0),c(0.0,4.71990005675,1.80806421444,0.0,0.0,1.59501894783,0.0,0.0),c(3.63905747464,5.39417149343,3.61612842889,4.24698990709,1.50824129844,2.39252842175,2.04870137171,0.0),c(4.5488218433,2.02281431004,2.71209632167,10.6174747677,0.0,0.0,2.04870137171,0.0))
targetgene="G0308"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(10.0074080553,10.1140715502,11.7524173939,12.7409697213,12.0659303875,11.9626421087,5.12175342928,8.04019759263),c(0.90976436866,1.34854287336,3.61612842889,10.6174747677,1.50824129844,0.797509473916,2.04870137171,2.68006586421),c(11.8269367926,8.09125724014,10.8483852867,8.49397981418,4.52472389531,11.9626421087,6.14610411513,5.36013172842),c(10.0074080553,8.09125724014,10.8483852867,8.49397981418,1.50824129844,6.38007579132,3.07305205757,8.04019759263),c(4.5488218433,2.69708574671,13.5604816083,4.24698990709,6.03296519374,3.19003789566,5.12175342928,0.0),c(5.45858621196,5.39417149343,9.94435317944,8.49397981418,6.03296519374,5.58256631741,3.07305205757,2.68006586421),c(4.5488218433,7.41698580346,4.52016053611,6.37048486064,9.04944779061,6.38007579132,10.2435068586,8.04019759263),c(4.5488218433,3.37135718339,6.32822475055,6.37048486064,3.01648259687,2.39252842175,2.04870137171,5.36013172842),c(2.72929310598,2.69708574671,4.52016053611,6.37048486064,6.03296519374,0.797509473916,1.02435068586,5.36013172842),c(4.5488218433,8.76552867682,9.04032107222,10.6174747677,3.01648259687,7.17758526524,3.07305205757,8.04019759263))
targetgene="G0620"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2.72929310598,3.37135718339,2.71209632167,2.12349495355,0.0,0.797509473916,0.0,5.36013172842),c(1.81952873732,3.37135718339,0.904032107222,0.0,1.50824129844,0.797509473916,0.0,0.0),c(21.8343448478,20.2281431004,9.04032107222,21.2349495355,9.04944779061,7.97509473916,11.2678575444,5.36013172842),c(0.90976436866,1.34854287336,1.80806421444,2.12349495355,1.50824129844,0.0,0.0,0.0),c(2.72929310598,4.71990005675,2.71209632167,4.24698990709,0.0,3.19003789566,0.0,5.36013172842),c(7.27811494928,3.37135718339,7.23225685778,2.12349495355,1.50824129844,1.59501894783,1.02435068586,0.0),c(1.81952873732,2.02281431004,2.71209632167,0.0,1.50824129844,1.59501894783,0.0,0.0),c(0.90976436866,4.04562862007,1.80806421444,2.12349495355,0.0,2.39252842175,1.02435068586,0.0),c(0.90976436866,0.674271436679,0.0,2.12349495355,3.01648259687,0.0,1.02435068586,0.0),c(0.0,2.69708574671,0.904032107222,0.0,0.0,0.0,2.04870137171,0.0))
targetgene="G0884"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(6.36835058062,8.09125724014,6.32822475055,10.6174747677,13.5741716859,9.57011368699,5.12175342928,2.68006586421),c(10.0074080553,1.34854287336,13.5604816083,6.37048486064,12.0659303875,3.98754736958,17.4139616595,2.68006586421),c(7.27811494928,5.39417149343,6.32822475055,16.9879596284,10.557689089,6.38007579132,7.17045480099,8.04019759263),c(6.36835058062,14.1597001702,14.4645137156,19.1114545819,13.5741716859,11.9626421087,12.2922082303,13.400329321),c(5.45858621196,2.69708574671,4.52016053611,6.37048486064,0.0,2.39252842175,2.04870137171,0.0),c(4.5488218433,4.04562862007,6.32822475055,6.37048486064,3.01648259687,4.78505684349,3.07305205757,2.68006586421),c(13.6464655299,11.4626144235,12.6564495011,10.6174747677,9.04944779061,13.5576610566,11.2678575444,13.400329321),c(5.45858621196,5.39417149343,9.94435317944,6.37048486064,0.0,0.797509473916,7.17045480099,5.36013172842),c(9.0976436866,2.69708574671,9.94435317944,6.37048486064,6.03296519374,2.39252842175,6.14610411513,0.0),c(3.63905747464,3.37135718339,7.23225685778,10.6174747677,4.52472389531,2.39252842175,12.2922082303,0.0))
targetgene="G0113"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(0.90976436866,4.71990005675,3.61612842889,4.24698990709,1.50824129844,2.39252842175,4.09740274342,0.0),c(0.0,4.04562862007,0.0,0.0,0.0,0.797509473916,1.02435068586,0.0),c(2.72929310598,6.06844293011,4.52016053611,0.0,0.0,2.39252842175,0.0,0.0),c(5.45858621196,6.74271436679,2.71209632167,2.12349495355,0.0,0.797509473916,2.04870137171,0.0),c(1.81952873732,2.02281431004,2.71209632167,0.0,4.52472389531,0.0,1.02435068586,0.0),c(0.90976436866,3.37135718339,4.52016053611,2.12349495355,1.50824129844,3.98754736958,2.04870137171,0.0),c(9.0976436866,8.09125724014,4.52016053611,2.12349495355,1.50824129844,7.97509473916,0.0,0.0),c(5.45858621196,2.69708574671,5.42419264333,6.37048486064,1.50824129844,0.797509473916,0.0,0.0),c(1.81952873732,2.02281431004,3.61612842889,10.6174747677,4.52472389531,3.19003789566,3.07305205757,0.0),c(0.0,2.02281431004,1.80806421444,2.12349495355,0.0,0.797509473916,0.0,2.68006586421))
targetgene="G0295"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=9
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("G0825","G0052","G0570","G0491","G0845","G0422","G0280","G0942","G0567","G0799")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='H1,H2,H3,H4_vs_L1,L2,L3,L4 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(24.5636379538,18.879600227,30.7370916455,33.9759192567,375.55208331,245.632917966,453.787353834,281.406915742),c(63.6835058062,51.2446291876,38.8733806105,63.7048486064,668.150895207,829.409852872,503.980537441,739.698178522),c(63.6835058062,44.5019148208,39.7774127178,59.4578586993,668.150895207,541.508932789,463.006510007,702.177256423),c(178.313816257,211.04695968,159.109650871,176.250081144,2328.72456478,2880.60421978,1991.3377333,1969.84841019),c(42.758925327,37.759200454,36.1612842889,44.5933940245,446.439424337,408.324850645,615.634762199,546.733436299),c(3.63905747464,10.7883429869,1.80806421444,6.37048486064,69.379099728,129.196534774,2.04870137171,182.244478766),c(43.6686896957,28.9936717772,33.4491879672,29.7289293496,561.065763018,359.676772736,356.474038678,458.29126278),c(47.3077471703,34.3878432706,33.4491879672,38.2229091638,588.21410639,403.539793801,367.741896222,428.810538273),c(26.3831666911,30.3422146505,22.6008026805,44.5933940245,342.370774745,327.776393779,254.038970092,627.135412225),c(26.3831666911,43.8276433841,28.0249953239,29.7289293496,472.07952641,557.459122267,457.884756577,420.770340681))
targetgene="G0825"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(26.3831666911,25.6223145938,24.408866895,44.5933940245,235.285642556,374.031943266,321.646115359,391.289616174),c(5.45858621196,5.39417149343,6.32822475055,10.6174747677,76.9203062202,73.3708716002,68.6314959523,45.5611196915),c(28.2026954284,27.6451289038,41.5854769322,91.3102830024,263.942227226,332.561450623,371.839298965,659.296202595),c(17.2855230045,6.74271436679,15.3685458228,8.49397981418,150.824129844,99.6886842395,187.456175512,37.5209220989),c(66.4127989122,62.7072436111,55.1459585405,72.1988284205,683.233308191,584.57444438,583.879890937,538.693238706),c(55.4956264882,60.0101578644,45.2016053611,31.8524243032,419.291080965,575.004330693,354.425337306,367.169023397),c(5.45858621196,8.09125724014,1.80806421444,8.49397981418,6.03296519374,4.78505684349,2.04870137171,8.04019759263),c(7.27811494928,11.4626144235,15.3685458228,8.49397981418,64.8543758327,80.5484568655,157.750005622,80.4019759263),c(9.0976436866,2.02281431004,12.6564495011,23.358444489,120.659303875,20.7352463218,102.435068586,251.926191236),c(10.9171724239,19.5538716637,10.8483852867,2.12349495355,107.085132189,128.3990253,120.873380931,50.92125142))
targetgene="G0052"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(25.4734023225,26.2965860305,18.9846742517,14.8644646748,200.596092692,197.782349531,123.946432988,123.283029754),c(34.5710460091,36.4106575806,20.7927384661,63.7048486064,321.255396567,282.318353766,295.012997526,348.408562347),c(4.5488218433,8.76552867682,4.52016053611,14.8644646748,46.7554802515,66.193286335,30.7305205757,29.4807245063),c(38.2101034837,24.2737717204,31.6411237528,31.8524243032,327.28836176,129.994044248,158.774356308,230.485664322),c(19.1050517419,12.8111572969,20.7927384661,14.8644646748,174.955990618,78.1559284437,167.99351248,168.844149445),c(7.27811494928,12.8111572969,2.71209632167,12.7409697213,67.8708584296,121.221440035,25.6087671464,67.0016466052),c(7.27811494928,8.76552867682,2.71209632167,4.24698990709,67.8708584296,64.5982673872,48.1444822352,26.8006586421),c(8.18787931794,15.5082430436,4.52016053611,10.6174747677,57.3131693405,103.676231609,43.0227288059,48.2411855558),c(4.5488218433,8.76552867682,3.61612842889,4.24698990709,19.6071368797,55.0281537002,47.1201315493,75.0418441978),c(8.18787931794,4.04562862007,9.04032107222,8.49397981418,57.3131693405,38.2804547479,88.0941589835,42.8810538273))
targetgene="G0570"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(42.758925327,46.5247291308,46.1056374683,55.2108687922,387.618013698,403.539793801,245.844164605,498.492250743),c(0.90976436866,6.74271436679,0.904032107222,4.24698990709,1.50824129844,91.7135895003,9.2191561727,67.0016466052),c(19.1050517419,16.1825144803,21.6967705733,38.2229091638,263.942227226,165.084461101,189.504876883,134.00329321),c(15.4659942672,17.5310573536,18.9846742517,23.358444489,120.659303875,180.237141105,172.090915224,131.323227346),c(20.9245804792,25.6223145938,13.5604816083,10.6174747677,294.107053195,204.959934796,109.605523387,115.242832161),c(3.63905747464,9.4398001135,6.32822475055,2.12349495355,16.5906542828,48.6480779089,19.4626630312,56.2813831484),c(30.9319885344,60.0101578644,35.2572521817,23.358444489,310.697707478,464.948023293,152.628252192,281.406915742),c(10.0074080553,13.4854287336,16.27257793,19.1114545819,90.4944779061,66.9907958089,115.751627502,85.7621076547),c(56.4053908569,34.3878432706,24.408866895,53.0873738386,429.848770054,324.586355884,257.11202215,418.090274817),c(25.4734023225,26.9708574671,26.2169311094,29.7289293496,144.79116465,188.212235844,221.259748145,278.726849878))
targetgene="G0491"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(15.4659942672,14.1597001702,13.5604816083,10.6174747677,102.560408294,52.6356252784,65.5584438947,88.4421735189),c(26.3831666911,28.9936717772,25.3128990022,19.1114545819,120.659303875,148.336762148,124.970783674,85.7621076547),c(27.2929310598,7.41698580346,10.8483852867,12.7409697213,101.052166995,52.6356252784,53.2662356645,58.9614490126),c(13.6464655299,16.1825144803,13.5604816083,16.9879596284,55.8049280421,79.7509473916,92.191561727,50.92125142),c(0.90976436866,3.37135718339,6.32822475055,8.49397981418,4.52472389531,7.97509473916,54.2905863503,45.5611196915),c(30.0222241658,22.2509574104,23.5048347878,8.49397981418,90.4944779061,65.3957768611,169.017863166,155.443820124),c(49.1272759076,39.1077433274,46.1056374683,27.6054343961,221.71147087,241.645370596,185.40747414,270.686652285),c(16.3757586359,8.09125724014,8.136288965,10.6174747677,48.2637215499,55.8256631741,32.7792219474,40.2009879631),c(51.8565690136,62.7072436111,49.7217658972,70.075333467,271.483433718,361.271791684,348.279233191,501.172316607),c(2.72929310598,7.41698580346,5.42419264333,6.37048486064,19.6071368797,35.8879263262,26.6331178322,37.5209220989))
targetgene="G0845"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(4.5488218433,12.1368858602,7.23225685778,6.37048486064,42.2307563562,67.7883052828,34.8279233191,53.6013172842),c(29.1124597971,23.5995002837,34.3532200744,67.9518385135,173.44774932,189.807254792,194.626630312,399.329813767),c(15.4659942672,15.5082430436,9.94435317944,10.6174747677,70.8873410265,111.651326348,81.9480548684,24.1205927779),c(32.7515172717,14.8339716069,22.6008026805,23.358444489,152.332371142,92.5110989742,135.214290533,155.443820124),c(30.9319885344,33.0393003972,30.7370916455,27.6054343961,173.44774932,231.277747436,187.456175512,93.8023052473),c(3.63905747464,2.69708574671,14.4645137156,12.7409697213,81.4450301155,24.7227936914,47.1201315493,75.0418441978),c(8.18787931794,14.8339716069,18.9846742517,8.49397981418,43.7389976546,104.473741083,89.1185096694,72.3617783336),c(22.7441092165,22.9252288471,19.8887063589,12.7409697213,140.266440754,145.146724253,122.922082303,99.1624369757),c(0.90976436866,0.0,0.904032107222,2.12349495355,4.52472389531,1.59501894783,2.04870137171,0.0),c(10.9171724239,6.74271436679,24.408866895,21.2349495355,85.9697540108,45.4580400132,214.089293344,128.643161482))
targetgene="G0422"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(27.2929310598,41.8048290741,32.54515586,27.6054343961,114.626338681,136.37412004,116.775978187,147.403622531),c(7.27811494928,0.0,12.6564495011,2.12349495355,45.2472389531,20.7352463218,43.0227288059,8.04019759263),c(27.2929310598,22.2509574104,27.1209632167,8.49397981418,155.348853739,101.283703187,143.40909602,107.202634568),c(43.6686896957,37.759200454,23.5048347878,27.6054343961,138.758199456,172.262046366,100.386367214,123.283029754),c(11.8269367926,8.76552867682,20.7927384661,12.7409697213,60.3296519374,51.8381158045,82.9724055543,67.0016466052),c(30.9319885344,37.0849290173,22.6008026805,29.7289293496,123.675786472,132.38657267,141.360394648,96.4823711115),c(22.7441092165,28.3194003405,9.94435317944,14.8644646748,90.4944779061,80.5484568655,53.2662356645,123.283029754),c(19.1050517419,16.856785917,24.408866895,12.7409697213,72.3955823249,95.7011368699,89.1185096694,37.5209220989),c(30.9319885344,8.09125724014,17.1766100372,23.358444489,132.725234262,39.0779642219,62.4853918372,152.76375426),c(5.45858621196,8.76552867682,6.32822475055,2.12349495355,10.557689089,61.4082294915,14.340909602,13.400329321))
targetgene="G0280"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(16.3757586359,6.06844293011,9.04032107222,10.6174747677,73.9038236233,36.6854358001,68.6314959523,42.8810538273),c(43.6686896957,26.9708574671,32.54515586,31.8524243032,137.249958158,174.654574788,154.676953564,163.484017717),c(94.6154943406,114.626144235,107.579820759,110.421737584,431.357011352,569.421764376,436.373392174,490.45205315),c(13.6464655299,18.2053287903,18.0806421444,16.9879596284,92.0027192045,103.676231609,83.9967562401,101.84250284),c(19.1050517419,26.2965860305,20.7927384661,23.358444489,98.0356843983,153.919328466,80.9237041826,134.00329321),c(3.63905747464,10.7883429869,12.6564495011,6.37048486064,1.50824129844,12.7601515826,11.2678575444,8.04019759263),c(12.7367011612,21.5766859737,13.5604816083,14.8644646748,70.8873410265,112.448835822,82.9724055543,72.3617783336),c(10.0074080553,10.1140715502,7.23225685778,25.4819394425,43.7389976546,55.0281537002,35.8522740049,72.3617783336),c(3.63905747464,0.0,8.136288965,2.12349495355,33.1813085656,2.39252842175,45.0714301776,18.7604610495),c(27.2929310598,33.7135718339,37.0653163961,21.2349495355,123.675786472,94.106117922,212.040591972,150.083688396))
targetgene="G0942"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(7.27811494928,2.69708574671,0.904032107222,4.24698990709,45.2472389531,9.57011368699,9.2191561727,8.04019759263),c(4.5488218433,0.674271436679,5.42419264333,4.24698990709,19.6071368797,6.38007579132,15.3652602878,8.04019759263),c(19.1050517419,21.5766859737,23.5048347878,31.8524243032,73.9038236233,79.7509473916,69.6558466382,56.2813831484),c(11.8269367926,16.1825144803,4.52016053611,25.4819394425,39.2142737593,63.0032484393,40.9740274342,40.2009879631),c(6.36835058062,10.7883429869,10.8483852867,4.24698990709,37.7060324609,33.4953979045,53.2662356645,16.0803951853),c(9.0976436866,16.856785917,26.2169311094,4.24698990709,31.6730672671,49.4455873828,32.7792219474,48.2411855558),c(13.6464655299,5.39417149343,12.6564495011,8.49397981418,31.6730672671,21.5327557957,58.3879890937,24.1205927779),c(12.7367011612,18.879600227,24.408866895,14.8644646748,36.1977911624,70.9783431785,51.2175342928,53.6013172842),c(23.6538735852,19.5538716637,18.0806421444,14.8644646748,110.101614786,68.5858147567,50.1931836069,34.8408562347),c(11.8269367926,18.879600227,21.6967705733,19.1114545819,31.6730672671,61.4082294915,63.509742523,61.6415148768))
targetgene="G0567"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(29.1124597971,20.902414537,23.5048347878,46.716888978,179.480714514,150.72929057,111.654224758,155.443820124),c(10.0074080553,10.1140715502,6.32822475055,4.24698990709,6.03296519374,13.5576610566,8.19480548684,5.36013172842),c(48.217511539,35.0621147073,28.9290274311,42.4698990709,232.269159959,152.324309518,193.602279627,174.204281174),c(10.0074080553,6.06844293011,18.0806421444,8.49397981418,40.7225150577,39.0779642219,107.556822015,101.84250284),c(46.3979828016,47.8732720042,41.5854769322,53.0873738386,227.744436064,206.554953744,160.823057679,206.365071544),c(13.6464655299,2.02281431004,4.52016053611,2.12349495355,75.4120649218,15.1526800044,24.5844164605,29.4807245063),c(9.0976436866,2.69708574671,6.32822475055,0.0,37.7060324609,12.7601515826,27.6574685181,16.0803951853),c(6.36835058062,5.39417149343,10.8483852867,6.37048486064,28.6565846703,44.6605305393,67.6071452664,13.400329321),c(4.5488218433,6.74271436679,11.7524173939,12.7409697213,37.7060324609,38.2804547479,61.4610411513,42.8810538273),c(38.2101034837,32.3650289606,34.3532200744,33.9759192567,123.675786472,164.286951627,160.823057679,171.524215309))
targetgene="G0799"
collabel=c("L1","L2","L3","L4","H1","H2","H3","H4")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



dev.off()
Sweave("GROWTH_10005_100_10_0.1_0.01_summary.Rnw");
library(tools);

texi2dvi("GROWTH_10005_100_10_0.1_0.01_summary.tex",pdf=TRUE);

