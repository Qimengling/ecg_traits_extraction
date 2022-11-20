source("fea_extract_functions_for_ukb.R")
fs<-500
ecg<-paste(rep(c("R","S","T","RS","ST","J_up",
                 "R-area","S-area","T-area","t_T","t_S","t_R",
                 "t_ST","t_QT"),12),
           c(rep("I",14),rep("II",14),rep("III",14),rep("aVR",14),
             rep("aVL",14),rep("aVF",14),rep("V1",14),rep("V2",14),
             rep("V3",14),rep("V4",14),rep("V5",14),rep("V6",14)),sep="_")
#3544314_2.rdata
#1001188_2.rdata
load("6024894_2.rdata")
data<-NULL
for (j in 1:12)
{
  data<-cbind(data,sub[[j]])
}

len_data<-10
data1<-data
for (j in 1:12)
{
  for (m in 1:10)
  {
    data1[((m-1)*fs+1):(m*fs),j]<-data[((m-1)*fs+1):(m*fs),j]-mean(data[((m-1)*fs+1):(m*fs),j])
  }
}

data2<-data1
for (run in 1:5)
{
  tt<-floor(dim(data2)[1]/len_data)
  range_r<-matrix(nrow = len_data,ncol = 12)
  temp<-matrix(0,nrow = len_data,ncol = 12)
  for (i in 1:12)
  {
    for (j in 1:len_data)
    {
      range_r[j,i]<-(max(data2[((j-1)*tt+1):(j*tt),i])-min(data2[((j-1)*tt+1):(j*tt),i]))
    }
  }
  for (i in 1:12)
  {
    temp[which(range_r[,i]>=(median(range_r[which(range_r[,i]!=0),i])*2) | range_r[,i]==0),i]<-1
  }
  
  Trange<-which(rowSums(temp)>0)
  if (length(Trange)!=0)
  {
    index<-NULL
    for (i in 1:length(Trange))
    {
      index<-c(index,c(((Trange[i]-1)*tt+1):(Trange[i]*tt)))
    }
    data2<-data2[-index,]
  }
  
  if (dim(data2)[1]/dim(data1)[1] <= 0.5)
  {
    break
  }
}

if (dim(data2)[1]/dim(data1)[1] > 0.5)
{
  TS_lead7<-find_Ts(7,data2,fs)
  
  if ((length(TS_lead7)/dim(data2)[1])*dim(data1)[1]>=5)
  {
    TS_lead1<-locate_other_lead_negative(TS_lead7,data2[,1],50)
  }
  if ((length(TS_lead7)/dim(data2)[1])*dim(data1)[1]<5)
  {
    TS_lead1<-find_Ts(1,data2,fs)
  }
  if ((length(TS_lead1)/dim(data2)[1])*dim(data1)[1]<5)
  {
    next
  }
  index1<-which(TS_lead1-230<=0)
  index2<-which(TS_lead1+230>dim(data2)[1])
  index<-union(index1,index2)
  if (length(index)>0)
  {
    TS_lead1<-TS_lead1[-index]
  }
  temp<-matrix(nrow = length(TS_lead1),ncol = 461)
  for (m in 1:length(TS_lead1))
  {
    temp[m,]<-data2[(TS_lead1[m]-230):(TS_lead1[m]+230),1]
  }
  com_wave<-matrix(nrow = 12,ncol = 461)
  com_wave[1,]<-colMeans(temp)
  for (j in 2:12)
  {
    if (j==4)
    {
      TS<-locate_other_lead_positive(TS_lead1,data2[,j],50)
    }
    if (j!=4)
    {
      TS<-locate_other_lead_negative(TS_lead1,data2[,j],50)
    }
    index1<-which(TS-230<=0)
    index2<-which(TS+230>dim(data2)[1])
    index<-union(index1,index2)
    if (length(index)>0)
    {
      TS<-TS[-index]
    }
    temp<-matrix(nrow = length(TS),ncol = 461)
    for (m in 1:length(TS))
    {
      temp[m,]<-data2[(TS[m]-230):(TS[m]+230),j]
    } 
    com_wave[j,]<-colMeans(temp)
  }
  combine_wave<-c(com_wave[1,],com_wave[2,],
                  com_wave[3,],com_wave[4,],
                  com_wave[5,],com_wave[6,],
                  com_wave[7,],com_wave[8,],
                  com_wave[9,],com_wave[10,],
                  com_wave[11,],com_wave[12,])
  if (is.na(combine_wave[1])==FALSE)
  {
    combine_wave_feature<-single_wave_fea_extract(combine_wave,fs)
  }
  
  names(combine_wave_feature)<-ecg
}
plot(data[,1],type = "l")
plot(data2[,7],type = "l")
plot(combine_wave,type = "l")

