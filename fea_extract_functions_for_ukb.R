single_wave_fea_extract<-function(data,fs)
{
result<-NULL
data1<-matrix(nrow=12,ncol=461)
for (i in 1:12)
{
data1[i,]<-data[((i-1)*461+1):(i*461)]
}

for (i in 1:12)
{
if (i==4)
{
S_pos<-c((231-50):(231+50))[which.max(data1[i,(231-50):(231+50)])]
R_pos<-c((S_pos-50):S_pos)[which.min(data1[i,(S_pos-50):S_pos])]
temp_j<--((data1[i,(S_pos+1):(S_pos+150)]-data1[i,S_pos:(S_pos+150-1)]))
j<-2
while (temp_j[j]>=mean(temp_j[1:(j-1)])/4 | temp_j[j+1]>=mean(temp_j[1:j])/4) 
{
  j<-j+1
  if ((j+1)>length(temp_j))
  {
    break
  }
}
J<-S_pos+j-1
temp_j<--(data1[i,R_pos:(R_pos-100+1)]-data1[i,(R_pos-1):(R_pos-100)])
j<-2
  while (temp_j[j]>=mean(temp_j[1:(j-1)])/4 | temp_j[j+1]>=mean(temp_j[1:j])/4) 
  {
    j<-j+1
    if (j+1>length(temp_j))
    {
      break
    }
  }
  R_start<-R_pos-j+1

}
if (i!=4)
{
S_pos<-c((231-50):(231+50))[which.min(data1[i,(231-50):(231+50)])]
R_pos<-c((S_pos-50):S_pos)[which.max(data1[i,(S_pos-50):S_pos])]
temp_j<-data1[i,(S_pos+1):(S_pos+150)]-data1[i,S_pos:(S_pos+150-1)]
j<-2
while (temp_j[j]>=mean(temp_j[1:(j-1)])/4 | temp_j[j+1]>=mean(temp_j[1:j])/4) 
{
  j<-j+1
  if ((j+1)>length(temp_j))
  {
    break
  }
}
J<-S_pos+j-1
temp_j<-data1[i,R_pos:(R_pos-100+1)]-data1[i,(R_pos-1):(R_pos-100)]
j<-2
  while (temp_j[j]>=mean(temp_j[1:(j-1)])/4 | temp_j[j+1]>=mean(temp_j[1:j])/4) 
  {
    j<-j+1
    if (j+1>length(temp_j))
    {
      break
    }
  }
  R_start<-R_pos-j+1
}
temp<-c(S_pos:S_pos+150)
if (max(data1[i,temp])>=abs(min(data1[i,temp])))
{
T_pos<-temp[which.max(data1[i,temp])]
}
if (max(data1[i,temp])<abs(min(data1[i,temp])))
{
T_pos<-temp[which.min(data1[i,temp])]
}

R<-data1[i,R_pos]
S<-abs(data1[i,S_pos])
T<-abs(data1[i,T_pos])
T_direction<-ifelse(1,-1,data1[i,T_pos]>0)

RS<-data1[i,R_pos]-data1[i,S_pos]
ST<-data1[i,T_pos]-data1[i,S_pos]
J_up<-data1[i,J]-data1[i,R_start]

d1<-data1[i,seq(R_pos,R_pos-100,by=-1)]
d2<-data1[i,seq(R_pos,S_pos,by=1)]

d3<-data1[i,seq(S_pos,R_pos,by=-1)]
d4<-data1[i,seq(S_pos,T_pos,by=1)]

d5<-data1[i,seq(T_pos,S_pos,by=-1)]
d6<-data1[i,seq(T_pos,461,by=1)]
if (i!=4)
{
TR1<-which(d1<=0)
TR2<-which(d2<=0)
TS1<-which(d3>=0)
TS2<-which(d4>=0)

    if (length(TR1)==0)
    {
      TR1<-which.min(d1)
    }
    
    if (length(TR2)==0)
    {
      TR2<-which.min(d2)
    }
    
    if (length(TS1)==0)
    {
      TS1<-which.max(d3)
    }
    
    if (length(TS2)==0)
    {
      TS2<-which.max(d4)
    }
}
if (i==4)
{
TR1<-which(d1>=0)
TR2<-which(d2>=0)
TS1<-which(d3<=0)
TS2<-which(d4<=0)

    if (length(TR1)==0)
    {
      TR1<-which.max(d1)
    }
    
    if (length(TR2)==0)
    {
      TR2<-which.max(d2)
    }
    
    if (length(TS1)==0)
    {
      TS1<-which.min(d3)
    }
    
    if (length(TS2)==0)
    {
      TS2<-which.min(d4)
    }
}

if (T_direction==-1)
{
  TA1<-which(d5>=0)
  TA2<-which(d6>=0)
  if (length(TA1)==0)
  {
  TA1<-which.max(d5)
  }
  if (length(TA2)==0)
  {
  TA2<-which.max(d6)
  }
}
if (T_direction==1)
{
  TA1<-which(d5<=0)
  TA2<-which(d6<=0)
  if (length(TA1)==0)
  {
  TA1<-which.min(d5)
  }
  if (length(TA2)==0)
  {
  TA2<-which.min(d6)
  }
}

RL<-R_pos-TR1[1]+1
RR<-R_pos+TR2[1]-1
SL<-S_pos-TS1[1]+1
SR<-S_pos+TS2[1]-1
AL<-T_pos-TA1[1]+1
AR<-T_pos+TA2[1]-1

inteR<-sum(data1[i,RL:RR])/fs
inteS<-sum(abs(data1[i,SL:SR]))/fs
inteA<-sum(abs(data1[i,AL:AR]))/fs

t_a<-(AR-AL)/fs
t_s<-(SR-SL)/fs
t_r<-(RR-RL)/fs

t_st<-(AR-SL)/fs
t_qt<-(AR-RL)/fs

result<-c(result,R,S,T,RS,ST,J_up,inteR,inteS,inteA,t_a,t_s,t_r,t_st,t_qt)
}
return(result)
}   
  
find_peak<-function(data,t,thres)
{
  index<-NULL
  for (i in 2:(length(data)-1))
  {
    if (data[i-1]<=data[i] & data[i]>=data[i+1] & data[i]> thres)
    {
      index<-c(index,i)
  } }
  return(t[index])
}

find_Ts<-function(lead_index,data3,fs)
{
  ####V1导联
  ###locate S wave
  x_recon<-data3[,lead_index]
  lendata<-length(x_recon)
  t<-1:lendata
  thres_s<-(min(x_recon[1:round(lendata*0.1)])+
              min(x_recon[round(lendata*0.1+1):round(lendata*0.2)])+
              min(x_recon[round(lendata*0.2+1):round(lendata*0.3)])+
              min(x_recon[round(lendata*0.3+1):round(lendata*0.4)])+
              min(x_recon[round(lendata*0.4+1):round(lendata*0.5)])+
              min(x_recon[round(lendata*0.5+1):round(lendata*0.6)])+
              min(x_recon[round(lendata*0.6+1):round(lendata*0.7)])+
              min(x_recon[round(lendata*0.7+1):round(lendata*0.8)])+
              min(x_recon[round(lendata*0.8+1):round(lendata*0.9)])+
              min(x_recon[round(lendata*0.9+1):lendata]))/30
  Ts<-find_peak(-x_recon,t,-thres_s)
  s<-x_recon[Ts]
  ###filter S wave
  if (length(Ts) > 100)
  {
    return(0)
  }
  Ts<-filter_by_value_and_range_ukb(Ts,x_recon,thres_s)
  s<-x_recon[Ts]
  
  ###
  TS_lead<-NULL
  for (i in 1:length(Ts))
  {
    if (Ts[i]-fs/2>0 & Ts[i]+fs/2<=lendata)
    {
      TS_lead<-c(TS_lead,Ts[i])
    }
  }
  return(TS_lead)
}

filter_by_value_and_range_ukb<-function(Tv,lead_data,thres)
{
  v<-lead_data[Tv]
  temp<-v[order(v,decreasing = T)]
  index3<-which(v>=mean(temp[3:length(temp)-2])/2)
  if (length(index3)>0)
  {
    Tv<-Tv[-index3]
  }
  v<-lead_data[Tv]
  
  for (run in 1:10)
  {
    A<-rep(NA,(length(Tv)-1))
    
    for (j in 2:length(Tv))
    {
      A[j-1]<-Tv[j]-Tv[j-1]
    }
    A1<-A
    if (length(which(A<20))>0)
    {
      A1<-A[-which(A<20)]
    }
    temp<-A1[order(A1,decreasing = T)]
    B<-median(temp[2:(length(A1)-1)])
    k<-NULL
    a<-NULL
    for (i in 1:length(A))
    {
      if (A[i]<(0.5*B))
      {
        if (v[i]>v[i+1])
        {
          k <-c(k,i)
        }
        if (v[i]<=v[i+1])
        {
          k <-c(k,(i+1))
        }
      }
      if (A[i]>(1.5*B))
      {
        temp<-(Tv[i]+1):(Tv[i+1]-1)
        temp1<-find_peak(-lead_data[temp],temp,0)
        if (length(temp1)>0)
        {
          a<-c(a,temp1[which.min(lead_data[temp1])])
        }
      }  
    }
    if (length(k)>0)
    {
      Tv<-Tv[-k]
    }
    Tv<-c(Tv,a)
    Tv<-Tv[order(Tv,decreasing = F)]
    v<-lead_data[Tv]
    temp<-which(v>thres)
    if (length(temp>0))
    {
      Tv<-Tv[-temp]
      v<-lead_data[Tv]
    }
  }
  return(Tv)
}

locate_other_lead_negative<-function(Tv,lead_data,flank_len)
{
  TV<-NULL
  for (i in 1:length(Tv))
  {
    if (Tv[i]-flank_len>0 & Tv[i]+flank_len<length(lead_data))
    {
      TV<-c(TV,Tv[i]-flank_len-1+which.min(lead_data[(Tv[i]-flank_len):(Tv[i]+flank_len)]))
     }
  }
  return(TV)
}

locate_other_lead_positive<-function(Tv,lead_data,flank_len)
{
  TV<-NULL
  for (i in 1:length(Tv))
  {
    if (Tv[i]-flank_len>0 & Tv[i]+flank_len<length(lead_data))
    {
      TV<-c(TV,Tv[i]-flank_len-1+which.max(lead_data[(Tv[i]-flank_len):(Tv[i]+flank_len)]))
    }
  }
  return(TV)
}