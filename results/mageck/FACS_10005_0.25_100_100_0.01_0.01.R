pdf(file='FACS_10005_0.25_100_100_0.01_0.01.pdf',width=4.5,height=4.5);
gstable=read.table('FACS_10005_0.25_100_100_0.01_0.01.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("G0991","G0293","G0099","G0407","G0520","G0380","G0689","G0135","G0904","G0810")
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
targetmat=list(c(494.185568819,737.288506082,666.081112994,504.935622399,0.0,8.7054158376,9.31293578597,20.4156005928),c(567.202280271,424.471480347,300.055709048,466.17288775,2.00652148358,11.6072211168,10.3477064289,5.56789107076),c(194.711230541,340.499948366,260.686697577,176.472449849,14.0456503851,0.9672684264,4.13908257154,0.0),c(362.967130845,352.495881506,258.558642903,253.997919146,27.0880400283,18.3781001016,10.3477064289,4.63990922563),c(713.235703177,742.825090608,655.440839623,525.337061688,14.0456503851,18.3781001016,24.8344954293,12.0637639866),c(400.004593176,386.638152751,445.827454224,466.17288775,13.0423896433,11.6072211168,11.3824770717,19.4876187476),c(322.755028886,215.926796525,485.196465695,379.466770773,17.0554326104,0.0,26.904036715,9.27981845126),c(245.505464595,454.922695242,266.006834263,303.981445404,12.0391289015,0.0,5.17385321443,12.0637639866),c(929.111197907,589.646252048,728.85872588,710.990159217,41.1336904134,25.1489790864,31.0431192866,22.271564283),c(326.987881724,238.995898717,195.781030017,322.342740764,10.0326074179,9.672684264,12.4172477146,17.6316550574))
targetgene="G0991"
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
targetmat=list(c(156.615555,238.07313463,321.336255789,100.98712448,23.0749970612,39.6580054824,42.4255963583,12.0637639866),c(257.145809899,323.890194787,264.942806926,341.724108088,0.0,9.672684264,0.0,25.0555098184),c(171.430539933,135.646320894,209.613385399,212.174968604,0.0,0.0,8.27816514309,0.0),c(196.82765696,370.95116326,262.814752252,181.572809671,0.0,4.836342132,0.0,13.9197276769),c(192.594804122,136.569084981,134.067444468,91.8064767999,0.0,0.0,1.03477064289,12.9917458318),c(240.214398548,295.284508068,170.244373928,280.519790222,4.01304296716,10.6399526904,0.0,12.9917458318),c(110.054173784,84.8942960695,122.363143761,369.266051128,0.0,8.7054158376,4.13908257154,10.2078002964),c(1829.65063916,1350.00386032,1457.71745176,1288.35089109,27.0880400283,24.18171066,30.0083486437,28.7674371989),c(645.510057772,419.857659909,604.367527445,620.203754381,8.02608593431,8.7054158376,14.4867890004,11.1357821415),c(1558.74805754,1426.5932796,1590.72086889,2139.09090944,37.1206474462,67.708789848,37.2517431439,69.5986383845))
targetgene="G0293"
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
targetmat=list(c(248.680104223,237.150370542,271.326970948,299.901157546,8.02608593431,0.0,0.0,13.9197276769),c(447.624187602,480.760089698,486.260493032,655.906273137,6.01956445074,20.3126369544,26.904036715,10.2078002964),c(63.4927925677,174.402412578,167.052291917,147.910434844,14.0456503851,25.1489790864,19.6606422148,16.7036732123),c(118.51987946,81.2032397186,183.012701973,99.9670525154,0.0,6.7708789848,0.0,13.9197276769),c(602.123316184,502.906427803,704.386097128,530.43742151,8.02608593431,11.6072211168,7.2433945002,2.78394553538),c(557.678361386,568.42267803,254.302533555,785.455412621,19.061954094,15.4762948224,25.8692660722,18.5596369025),c(456.089893278,450.308874803,738.434971913,520.236701866,10.0326074179,33.854394924,16.5563302862,0.0),c(153.440915372,356.186937857,284.095298993,303.981445404,8.02608593431,9.672684264,0.0,0.0),c(191.536590913,230.691021928,34.0488747856,226.455976106,11.0358681597,0.0,2.06954128577,3.7119273805),c(162.964834257,228.845493753,178.756592625,135.669571271,1.00326074179,0.0,9.31293578597,4.63990922563))
targetgene="G0099"
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
targetmat=list(c(279.368287298,262.987764998,337.296665845,144.850218951,3.00978222537,13.5417579696,4.13908257154,0.0),c(188.361951284,501.983663715,486.260493032,307.041661297,14.0456503851,5.8036105584,33.1126605724,4.63990922563),c(144.975209696,198.394278858,243.662260185,367.225907199,0.0,25.1489790864,12.4172477146,6.49587291588),c(590.48297088,655.162502275,554.358242604,566.139940266,4.01304296716,18.3781001016,28.9735780008,20.4156005928),c(193.653017331,61.8251938767,114.914952402,164.231586275,0.0,3.8690737056,8.27816514309,0.0),c(96.297402061,73.8211270169,127.683280446,59.1641739377,7.02282519253,0.0,3.10431192866,0.0),c(238.097972129,370.028399172,366.025403946,216.255256462,14.0456503851,21.2799053808,37.2517431439,11.1357821415),c(150.266275744,131.032500455,110.658843053,122.408635733,5.01630370895,0.0,0.0,0.0),c(691.013225778,609.947061977,891.654908449,886.442537101,9.0293466761,14.509026396,22.7649541435,23.1995461282),c(616.938301116,804.650284485,858.670061,676.307712426,22.0717363194,27.0835159392,20.6954128577,10.2078002964))
targetgene="G0407"
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
targetmat=list(c(321.696815676,393.097501365,523.501449829,413.129145599,6.01956445074,4.836342132,19.6606422148,19.4876187476),c(80.4242039191,91.3536446835,112.786897727,79.5656132265,7.02282519253,4.836342132,6.20862385732,1.85596369025),c(227.515840034,276.829226314,319.208201115,315.202237013,5.01630370895,11.6072211168,8.27816514309,8.35183660614),c(322.755028886,268.524349524,402.202333405,370.286123093,20.0652148358,19.345368528,22.7649541435,12.9917458318),c(233.865119291,316.508082085,235.150041488,228.496120035,14.0456503851,16.4435632488,0.0,19.4876187476),c(1110.06565673,1113.77625387,1036.36262629,1222.0462134,36.1173867044,66.7415214216,37.2517431439,27.8394553538),c(112.170600203,299.898328506,196.845057354,149.950578773,24.0782578029,0.0,12.4172477146,0.0),c(167.197687095,145.796725858,204.293248714,117.308275911,3.00978222537,16.4435632488,0.0,4.63990922563),c(242.330824967,338.65442019,240.470178174,318.262452906,0.0,5.8036105584,7.2433945002,24.1275279733),c(642.335418143,428.162536698,639.480429568,606.942818844,13.0423896433,29.9853212184,7.2433945002,25.9834916635))
targetgene="G0520"
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
targetmat=list(c(209.526215473,211.312976086,154.283963872,160.151298418,4.01304296716,0.0,5.17385321443,12.0637639866),c(169.314113514,23.0691021928,131.939389794,188.713313422,0.0,3.8690737056,0.0,4.63990922563),c(148.149849325,181.784525279,95.7624603346,142.810075022,5.01630370895,0.0,0.0,0.0),c(108.995960575,77.5121833678,18.0884647299,81.6057571554,0.0,0.0,0.0,5.56789107076),c(807.416678819,1147.91852511,1067.21941906,755.873325652,39.1271689298,42.5598107616,44.4951376441,46.3990922563),c(256.08759669,310.971497559,435.187180854,220.33554432,12.0391289015,0.0,2.06954128577,10.2078002964),c(382.014968616,403.24790633,661.825003646,630.404474026,6.01956445074,22.2471738072,6.20862385732,11.1357821415),c(169.314113514,210.390211998,188.332838658,257.05813504,5.01630370895,5.8036105584,5.17385321443,0.0),c(138.625930439,230.691021928,189.396865995,171.372090026,11.0358681597,0.0,0.0,4.63990922563),c(179.896245608,142.105669508,285.15932633,193.813673244,10.0326074179,11.6072211168,10.3477064289,12.9917458318))
targetgene="G0380"
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
targetmat=list(c(471.96309142,610.869826065,844.837705618,886.442537101,16.0521718686,21.2799053808,6.20862385732,29.695419044),c(1357.68754774,1352.77215259,1294.92126919,1303.65197056,72.2347734088,22.2471738072,73.4687156449,38.9752374953),c(233.865119291,197.47151477,378.79373199,371.306195057,3.00978222537,2.9018052792,4.13908257154,15.7756913671),c(196.82765696,165.1747717,185.140756647,207.074608782,15.0489111268,16.4435632488,0.0,0.0),c(200.002296588,331.272307489,126.619253109,292.760653795,0.0,24.18171066,0.0,11.1357821415),c(586.250118042,458.613751593,540.525887222,509.015910257,20.0652148358,19.345368528,1.03477064289,23.1995461282),c(722.759622062,1103.6258489,825.685213552,670.187280639,15.0489111268,18.3781001016,15.5215596433,26.9114735087),c(190.478377703,121.804859578,127.683280446,222.375688249,0.0,10.6399526904,13.4520183575,4.63990922563),c(33.8628227028,84.8942960695,109.594815716,0.0,5.01630370895,0.0,0.0,0.0),c(125.927371926,58.1341375258,75.5459409306,18.36129536,10.0326074179,0.0,0.0,0.0))
targetgene="G0689"
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
targetmat=list(c(14.8149849325,39.6788557716,38.3049841338,40.8028785777,26.0847792865,47.3961528936,41.3908257154,47.3270741014),c(178.838032399,165.1747717,197.909084691,194.833745209,126.410853465,169.27197462,186.258715719,181.884441645),c(39.1538887501,35.9877994208,98.9545423457,62.224389831,43.1402118969,26.1162475128,76.5730275736,35.2633101148),c(116.403453041,134.723556806,160.668127895,118.328347875,117.381506789,119.941284874,153.146055147,123.421585402),c(480.428797096,447.54058254,467.108000965,495.754974719,468.522766416,465.256113099,448.05568837,465.846886253),c(80.4242039191,99.6585214729,91.5063509864,123.428707698,104.339117146,95.7595742136,94.1641285026,127.133512782),c(34.9210359122,57.2113734381,44.6891481561,60.1842459021,31.1010829955,64.8069845688,66.2253211447,66.8146928491),c(66.6674321961,23.9918662805,13.8323553817,34.6824467911,65.2119482163,23.2144422336,4.13908257154,23.1995461282),c(121.694519088,50.7520248241,56.3934488637,75.4853253688,104.339117146,47.3961528936,48.6342202156,45.4711104112),c(0.0,31.3739789822,43.6251208191,20.4014392889,0.0,46.4288844672,20.6954128577,26.9114735087))
targetgene="G0135"
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
targetmat=list(c(206.351575845,179.016233016,179.820619962,89.766332871,195.635844649,181.846464163,167.632844148,116.925712486),c(31.7463962838,41.524383947,35.1129021227,72.4251094754,15.0489111268,38.690737056,20.6954128577,68.6706565393),c(91.0063360137,84.8942960695,143.643690502,92.8265487643,87.2836845357,88.0214268024,117.963853289,93.7261663577),c(83.5988435475,55.3658452627,74.4819135936,112.207916089,50.1630370895,50.2979581728,72.433945002,73.310565765),c(43.3867415879,46.1382043856,88.3142689752,35.7025187555,47.1532548641,42.5598107616,80.7121101451,31.5513827343),c(101.588468108,88.5853524203,102.146624357,67.3247496532,77.2510771178,66.7415214216,92.0945872169,87.2302934419),c(370.374623312,430.008064874,316.016119104,479.433823288,370.20321372,365.627465179,330.091835081,449.143213041),c(110.054173784,111.654454613,117.043007076,65.2846057243,135.440200142,93.8250373608,118.998623932,70.5266202296),c(29.6299698649,60.902429789,70.2258042454,56.1039580444,33.107604479,69.6433267008,48.6342202156,62.1747836235),c(111.112386993,125.495915929,73.4178862565,64.2645337599,87.2836845357,94.7923057872,66.2253211447,54.7509288624))
targetgene="G0904"
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
targetmat=list(c(15.8731981419,11.0731690525,19.1524920669,14.2810075022,10.0326074179,4.836342132,12.4172477146,14.847709522),c(142.858783277,92.2764087712,105.338706368,116.288203946,110.358681597,100.595916346,92.0945872169,116.925712486),c(170.372326723,319.276374348,181.948674636,177.492521813,193.629323165,294.049601626,190.397798291,225.499588366),c(491.01092919,353.418645594,413.906634113,449.851736319,490.594502735,311.460433301,384.934679154,425.015685068),c(38.0956755406,31.3739789822,58.5215035378,119.34841984,38.123908188,32.8871264976,66.2253211447,125.277549092),c(64.5510057772,62.7479579644,60.6495582119,67.3247496532,73.2380341506,59.0033740104,65.1905505018,78.8784568357),c(231.748692872,135.646320894,214.933522084,207.074608782,218.71084171,99.6286479193,167.632844148,143.837185995),c(39.1538887501,16.6097535788,64.9056675601,39.7828066133,27.0880400283,22.2471738072,65.1905505018,38.0472556502),c(207.409789054,161.48371535,272.390998285,201.97424896,151.49237201,151.861142945,224.545229506,133.629385698),c(64.5510057772,75.6666551924,91.5063509864,117.308275911,60.1956445074,99.6286479193,76.5730275736,93.7261663577))
targetgene="G0810"
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
targetgenelist=c("G0872","G0170","G0907","G0701","G0836","G0888","G0022","G0205","G0341","G0492")
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
targetmat=list(c(7.40749246623,1.84552817542,11.7043007076,4.08028785777,237.772795804,194.420953706,153.146055147,158.684895517),c(0.0,4.61382043856,0.0,1.02007196444,165.538022395,121.875821726,175.911009291,154.972968136),c(69.8420718245,35.065035333,40.4330388079,36.7225907199,1213.94549757,1423.81912366,999.588441028,1175.75299777),c(26.4553302365,1.84552817542,10.6402733705,0.0,354.151041852,200.224564265,320.778899295,282.106480918),c(0.0,0.0,0.0,0.0,186.606497973,117.039479594,219.371376292,103.005984809),c(23.2806906082,48.9064966487,52.1373395155,55.0838860799,992.22487363,960.497547416,1265.52449625,869.518988883),c(12.6985585135,3.69105635085,2.1280546741,0.0,154.502154236,169.27197462,152.111284504,115.997730641),c(7.40749246623,7.38211270169,4.2561093482,1.02007196444,349.134738143,423.663570763,259.727431364,240.347297888),c(0.0,6.45934861398,7.44819135936,0.0,80.2608593431,123.810358579,112.790000075,131.773422008),c(8.46570567569,4.61382043856,9.57624603346,9.18064767999,234.763013579,346.282096651,249.379724936,270.970698777))
targetgene="G0872"
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
targetmat=list(c(28.5717566555,16.6097535788,12.7683280446,15.3010794666,698.269476285,744.796688328,536.011193015,612.468017783),c(9.52391888515,25.8373944559,10.6402733705,13.2609355378,400.301035974,310.493164875,404.595321368,269.114735087),c(7.40749246623,0.0,0.0,15.3010794666,118.384767531,123.810358579,41.3908257154,235.707388662),c(0.0,0.0,0.0,0.0,138.449982367,144.122995534,49.6689908585,180.9564598),c(11.6403453041,7.38211270169,5.32013668525,20.4014392889,55.1793407984,47.3961528936,49.6689908585,100.222039274),c(2.11642641892,6.45934861398,10.6402733705,7.1405037511,192.626062424,320.165849139,417.012569083,226.427570211),c(10.5821320946,2.76829226314,0.0,0.0,336.092348499,301.787749037,140.728807432,258.90693479),c(2.11642641892,14.7642254034,0.0,2.04014392889,291.948875861,276.638769951,209.023669863,167.964713968),c(0.0,0.0,0.0,8.16057571554,187.609758715,237.948032895,190.397798291,163.324804742),c(0.0,6.45934861398,0.0,0.0,65.2119482163,166.370169341,191.432568934,270.970698777))
targetgene="G0170"
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
targetmat=list(c(17.9896245608,26.7601585436,11.7043007076,36.7225907199,937.045532831,784.454693811,1118.58706496,1013.35617488),c(12.6985585135,17.5325176665,0.0,3.06021589333,226.736927644,344.347559799,159.354679004,194.876187476),c(0.0,0.0,9.57624603346,20.4014392889,187.609758715,289.213259494,283.527156151,189.308296406),c(0.0,0.0,2.1280546741,0.0,86.2804237939,163.468364062,152.111284504,153.117004446),c(7.40749246623,4.61382043856,0.0,8.16057571554,258.841271382,218.602664367,166.598073505,149.405077065),c(0.0,3.69105635085,3.19208201115,9.18064767999,242.789099513,84.1523530968,169.702385433,224.571606521),c(102.646681318,108.88616235,130.875362457,114.248060018,233.759752837,214.733590661,323.883211223,295.09822675),c(0.0,2.76829226314,23.4086014151,4.08028785777,274.89344325,323.067654418,251.449266221,175.388568729),c(6.34927925677,0.0,0.0,23.4616551822,177.577151297,139.286653402,133.485412932,141.981222304),c(0.0,0.0,5.32013668525,0.0,146.476068301,89.9559636552,111.755229432,202.300042237))
targetgene="G0907"
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
targetmat=list(c(73.0167114528,79.3577115432,89.3782963123,56.1039580444,67.2184696999,62.872447716,90.0250459311,70.5266202296),c(22.2224773987,38.7560916839,30.8567927745,45.9032383999,24.0782578029,47.3961528936,27.9388073579,24.1275279733),c(31.7463962838,24.9146303682,4.2561093482,22.4415832177,32.1043437373,15.4762948224,10.3477064289,37.119273805),c(151.324488953,148.565018122,130.875362457,91.8064767999,158.515197203,163.468364062,122.102935861,100.222039274),c(30.6881830744,38.7560916839,29.7927654374,22.4415832177,52.169558573,46.4288844672,53.8080734301,22.271564283),c(58.2017265204,47.0609684733,51.0733121784,65.2846057243,44.1434726387,47.3961528936,28.9735780008,54.7509288624),c(125.927371926,117.191039139,127.683280446,123.428707698,147.479329043,93.8250373608,116.929082646,89.0862571321),c(82.540630338,107.963398262,74.4819135936,82.6258291199,93.3032489864,99.6286479193,65.1905505018,108.57387588),c(38.0956755406,78.4349474555,54.2653941896,19.3813673244,64.2086874745,76.4142056856,58.9819266445,33.4073464245),c(47.6195944258,27.6829226314,34.0488747856,56.1039580444,64.2086874745,35.7889317768,45.529908287,73.310565765))
targetgene="G0701"
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
targetmat=list(c(128.043798345,80.2804756309,82.99413229,105.067412338,108.352160113,86.0868899496,62.0862385732,97.4380937382),c(92.0645492232,34.1422712453,79.8020502788,71.405037511,94.3065097282,56.1015687312,78.6425688593,46.3990922563),c(146.033422906,158.715423086,139.387581154,146.89036288,159.518457945,155.73021665,139.69403679,160.540859207),c(85.7152699664,46.1382043856,60.6495582119,87.7261889421,68.2217304417,55.1343003048,61.0514679303,103.005984809),c(382.014968616,356.186937857,287.287381004,367.225907199,354.151041852,389.809175839,347.68293601,373.976683586),c(160.848407838,155.024366736,222.381713444,230.536263964,124.404331982,138.319384975,227.649541435,204.156005928),c(75.1331378718,128.264208192,118.107034413,107.107556267,88.2869452775,164.435632488,128.311559718,145.693149685),c(174.605179561,143.951197683,131.939389794,113.227988053,175.570629813,151.861142945,108.650917503,131.773422008),c(32.8046094933,26.7601585436,47.8812301673,43.863094471,52.169558573,46.4288844672,54.842844073,46.3990922563),c(42.3285283785,69.2073065784,44.6891481561,88.7462609065,55.1793407984,66.7415214216,84.8511927167,95.582130048))
targetgene="G0836"
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
targetmat=list(c(126.985585135,129.18697228,113.850925064,47.9433823288,147.479329043,147.992069239,117.963853289,72.3825839198),c(37.0374623312,34.1422712453,35.1129021227,16.3211514311,29.0945615119,33.854394924,26.904036715,15.7756913671),c(50.7942340542,32.2967430699,48.9452575043,71.405037511,70.2282519253,42.5598107616,65.1905505018,81.6624023711),c(183.070885237,131.032500455,138.323553817,149.950578773,175.570629813,115.104942742,115.894312003,144.76516784),c(39.1538887501,29.5284508068,56.3934488637,18.36129536,42.1369511552,57.0688371576,55.8776147158,34.3353282697),c(265.611515575,240.841426893,255.366560892,301.941301475,227.740188386,300.820480611,232.823394649,259.834916635),c(116.403453041,163.329243525,148.963827187,150.970650738,141.459764592,123.810358579,131.415871647,131.773422008),c(73.0167114528,28.6056867191,50.0092848414,45.9032383999,66.2152089581,54.1670318784,44.4951376441,36.1912919599),c(107.937747365,70.1300706661,112.786897727,116.288203946,112.36520308,52.2324950256,121.068165218,126.205530937),c(93.1227624326,100.581285561,105.338706368,116.288203946,79.2575986014,104.464990051,85.8859633595,107.645894035))
targetgene="G0888"
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
targetmat=list(c(159.790194629,112.577218701,134.067444468,138.729787164,152.495632752,122.843090153,156.250367076,154.972968136),c(71.9584982434,141.18290542,150.027854524,112.207916089,73.2380341506,184.748269442,155.215596433,116.925712486),c(101.588468108,85.8170601572,86.1862143011,72.4251094754,124.404331982,131.54850599,84.8511927167,81.6624023711),c(202.118723007,191.934930244,215.997549421,224.415832177,222.723884677,234.078959189,251.449266221,238.491334197),c(92.0645492232,79.3577115432,63.8416402231,52.0236701866,89.2902060193,65.7742529952,42.4255963583,41.7591830307),c(59.2599397298,57.2113734381,86.1862143011,39.7828066133,49.1597763477,69.6433267008,94.1641285026,50.1110196368),c(32.8046094933,30.4512148945,31.9208201115,24.4817271466,42.1369511552,41.5925423352,51.7385321443,37.119273805),c(92.0645492232,86.7398242449,82.99413229,51.0035982221,73.2380341506,92.8577689344,100.37275236,59.3908380881),c(155.557341791,186.398345718,139.387581154,122.408635733,123.40107124,190.551880001,143.833119361,124.349567247),c(76.1913510812,111.654454613,56.3934488637,95.8867646576,74.2412948924,141.221190254,62.0862385732,90.9422208224))
targetgene="G0022"
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
targetmat=list(c(20.1060509798,6.45934861398,12.7683280446,11.2207916089,23.0749970612,1.9345368528,13.4520183575,19.4876187476),c(22.2224773987,32.2967430699,42.561093482,56.1039580444,15.0489111268,46.4288844672,36.216972501,48.2550559466),c(61.3763661488,62.7479579644,10.6402733705,52.0236701866,71.231512667,76.4142056856,27.9388073579,54.7509288624),c(6.34927925677,5.53658452627,0.0,5.10035982221,8.02608593431,10.6399526904,0.0,0.927981845126),c(12.6985585135,95.967465122,39.3690114709,68.3448216177,14.0456503851,88.9886952288,34.1474312152,74.2385476101),c(101.588468108,146.719489946,123.427171098,76.5053973332,90.293466761,135.417579696,77.6077982165,81.6624023711),c(21.1642641892,36.9105635085,58.5215035378,24.4817271466,16.0521718686,29.018052792,67.2600917876,35.2633101148),c(32.8046094933,25.8373944559,39.3690114709,41.8229505422,22.0717363194,26.1162475128,51.7385321443,31.5513827343),c(69.8420718245,35.065035333,47.8812301673,34.6824467911,85.2771630521,61.9051792896,70.3644037162,34.3353282697),c(40.2121019595,13.8414613157,47.8812301673,17.3412233955,30.0978222537,26.1162475128,46.5646789299,36.1912919599))
targetgene="G0205"
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
targetmat=list(c(357.676064798,368.182870997,413.906634113,423.329865244,345.121695176,333.707607108,405.630092011,422.231739532),c(31.7463962838,16.6097535788,24.4726287522,41.8229505422,32.1043437373,16.4435632488,21.7301835006,42.6871648758),c(73.0167114528,49.8292607364,50.0092848414,72.4251094754,107.348899371,70.6105951272,57.9471560016,65.886711004),c(110.054173784,117.191039139,67.0337222342,102.007196444,79.2575986014,121.875821726,76.5730275736,129.917458318),c(119.578092669,59.0569016135,62.777612886,105.067412338,125.407592724,68.6760582744,67.2600917876,90.0142389772),c(157.67376821,133.800792718,188.332838658,121.388563769,147.479329043,133.483042843,175.911009291,140.125258614),c(34.9210359122,12.918697228,10.6402733705,56.1039580444,43.1402118969,16.4435632488,22.7649541435,57.5348743978),c(283.601140136,247.300775507,331.97652916,367.225907199,313.017351438,312.427701727,350.787247938,386.040447572),c(29.6299698649,58.1341375258,51.0733121784,27.54194304,20.0652148358,30.9525896448,50.7037615014,18.5596369025),c(47.6195944258,11.0731690525,32.9848474486,11.2207916089,44.1434726387,8.7054158376,42.4255963583,19.4876187476))
targetgene="G0341"
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
targetmat=list(c(252.912957061,304.512148945,260.686697577,227.476048071,274.89344325,236.980764468,243.171101078,227.355552056),c(31.7463962838,13.8414613157,25.5366560892,27.54194304,57.185862282,22.2471738072,48.6342202156,29.695419044),c(158.731981419,96.8902292097,113.850925064,183.6129536,185.603237231,105.432258478,140.728807432,178.172514264),c(197.885870169,107.963398262,85.1221869641,97.9269085865,151.49237201,155.73021665,99.3379817171,95.582130048),c(144.975209696,151.333310385,94.6984329975,127.508995555,164.534761653,143.155727107,101.407523003,148.47709522),c(64.5510057772,76.5894192801,50.0092848414,109.147700195,57.185862282,65.7742529952,65.1905505018,91.8702026675),c(148.149849325,152.256074472,74.4819135936,90.7864048354,167.544543879,120.9085533,91.059816574,82.5903842162),c(95.2391888515,105.195105999,139.387581154,122.408635733,107.348899371,117.039479594,142.798348718,108.57387588),c(242.330824967,235.304842366,215.997549421,202.994320924,227.740188386,231.17715391,180.050091862,239.419316043),c(40.2121019595,17.5325176665,29.7927654374,55.0838860799,44.1434726387,45.4616160408,42.4255963583,49.1830377917))
targetgene="G0492"
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
Sweave("FACS_10005_0.25_100_100_0.01_0.01_summary.Rnw");
library(tools);

texi2dvi("FACS_10005_0.25_100_100_0.01_0.01_summary.tex",pdf=TRUE);

