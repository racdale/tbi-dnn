library(ggplot2)
library(lme4)
library(stringdist)
library(stargazer)
library(gridExtra)

source('cogsci_functions.R')

a = read.csv('processed.csv',header=F)
colnames(a) = list('r','fl','i','j','who_i','who_j','n_i','n_j',
                   'txt_i','txt_j', 'gra_i','gra_j','h_1','h_2','baseline','layer')
# nb: raw text removed from data shared so as not to disseminate coelho corpus directly
# https://tbi.talkbank.org/access/English/Coelho.html
# nb: password protected, permission from talkbank.org easy to request
a[1,]

dim(a)

# assess transript size
docs = aggregate(i>0~fl,data=a,sum)
docs$`i > 0` = docs$`i > 0`/20/12
docs # cross ref to confirm
mean(docs$`i > 0`)

a$abs_dist = abs(a$i-a$j) # absolute distance across turns
a$rel_dist = a$i-a$j # relative distance
a = a[a$n_i>5&a$n_j>5,] # must be 6 or more tokens
a = a[a$abs_dist>0 & a$abs_dist<=10,] # no more than 10 away, but not the same turn

dim(a) # size post filtering

# assess transcript size post filtering
docs = aggregate(i>0~fl,data=a,sum)
docs$`i > 0` = docs$`i > 0`/20/12
docs # cross ref to confirm
mean(docs$`i > 0`)

# check sample to ensure most utterances clean (annotation, xxx, etc.)
a$txt_j[sample(1:nrow(a),100)]

a$lev = stringsim(a$txt_i,a$txt_j) # levenshtein distance, lexical
a$self = a$who_i==a$who_j # self = TRUE if same person
a$tb = FALSE 
a$tb[grep('/tb',a$fl)] = TRUE # tb = TRUE when TBI conversation

# factor by length of relevant utterance
a$h_1_n = as.numeric(a$h_1)/as.numeric(a$n_i) # focal variable (i)
a$h_2_n = as.numeric(a$h_2)/as.numeric(a$n_j)

a$residH = resid(lm(h_1_n~lev,data=a))
a$layerCat = as.factor(a$layer)
a$distC = a$abs_dist - mean(a$abs_dist)
a$par = (a$who_i =="*PAR") # participant = TRUE when i = TBI patient or control

sum(a$par[a$tb])
sum(a$par[!a$tb])

hist(a[a$par & a$layer==12,]$residH)

# let's use PAR only; the TB participant vs. control
processedData = a[a$par,]

# convergence issues; simplest random effect chosen
modl_base = lmer(residH~distC+(1|fl),data=processedData)
modl_tbdist = lmer(residH~tb*distC+(1|fl),data=processedData)
modl_tbselfdist = lmer(residH~self*tb*distC+(1|fl),data=processedData)
modl_tbselfdist_layer = lmer(residH~self*tb*layerCat*distC+(1|fl),data=processedData)

# check contributions of models compared to relative baselines
anova(modl_base,modl_tbdist)
anova(modl_base,modl_tbselfdist)
anova(modl_base,modl_tbselfdist_layer)

# get p's for chosen model
coefs = data.frame(coef(summary(modl_tbselfdist_layer)))
coefs$p = round(1-pnorm(abs(coefs$t.value)),5)
coefs
stargazer(coefs[,c(1,3:4)],digits=2,summary=F,rownames = T)

cor(predict(modl_base),processedData$residH)^2
cor(predict(modl_tbselfdist),processedData$residH)^2
cor(predict(modl_tbself_layer),processedData$residH)^2

# let's look at results across layers; build model for each
layer_results = c()
for (i in 1:12) {
  print(i) 
  
  modl_base_1L = lmer(residH~distC+(1|fl),data=processedData[processedData$layer==i,])
  modl_full = lmer(residH~self*tb*distC+(1|fl),data=processedData[processedData$layer==i,])
  
  rbase = cor(predict(modl_base_1L),processedData[processedData$layer==i,]$residH)^2
  rtbi = cor(predict(modl_full),processedData[processedData$layer==i,]$residH)^2
  pval = anova(modl_base_1L,modl_full)$P[2]
  chimod = anova(modl_base_1L,modl_full)$Ch[2]
  
  layer_results = rbind(layer_results,data.frame(rbase=rbase,rtbi=rtbi,
                                                 chimod=chimod,
                                                 layer=i,p=pval))
}

layer_results$rDiff = (layer_results$rtbi-layer_results$rbase)
stargazer(subset(layer_results,select=c(layer,chimod,rDiff,p)),
          digits=4,summary=F,rownames = F)
plot(layer_results$layer,layer_results$rDiff,type='o',lwd=2,xlab='BERT layer',ylab='Change in r-squared fit')

# plot temporal profiles over k and layer
targetLayer = 10
profile(1:4,T)
profile(1:4,F)

profile(9:12,T)
profile(9:12,F)

# first, build model for each individual
tbs = unique(processedData$fl[processedData$tb])
layerRes = c()
for (layer in 1:12) {
  print(layer)
  for (fl in tbs) {
    tmp = processedData[(!processedData$tb | processedData$fl==fl)&
                          processedData$layer==layer&
                          processedData$self&
                          processedData$abs_dist<5,]
    modl_1L = lmer(residH~tb*distC+(1|fl),
                   data=tmp)
    tv = coef(summary(modl_1L))[10]
    p = round(1-pnorm(abs(tv)),5)
    layerRes = rbind(layerRes, data.frame(layer=layer,
                                          fl=fl,
                                          t=tv,
                                          p=p))
  }
}

# let's store ts, ps, etc. for plotting
vecs = matrix(0,nrow=length(unique(layerRes$fl)),ncol=12)
for (layer in 1:12) {
  print(layer)
  for (fl in tbs) {
    p = layerRes[layerRes$layer==layer&layerRes$fl==fl,]$p
    participant = which(tbs==fl)
    s = sign(layerRes[layerRes$layer==layer&layerRes$fl==fl,]$t)
    vecs[participant,layer] = layerRes[layerRes$layer==layer&layerRes$fl==fl,]$t
  }
}

# order approximately by positivity
sums = rowSums(vecs*(vecs>0)-abs(vecs)*(vecs<0))
ixes = order(sums)

# plot individual differences
plot(c(1,nrow(vecs)),c(1,ncol(vecs)),cex=4,col='white',
     xlab='TBI participant (approximately ranked by coefficients)',ylab='BERT layer (L)',xaxt='n')
axis(side=1,at=c(1,24,length(tbs)))

for (layer in 1:12) {
  print(layer)
  for (j in ixes) {
    fl = tbs[j]
    t = layerRes[layerRes$layer==layer&layerRes$fl==fl,]$t
    p = layerRes[layerRes$layer==layer&layerRes$fl==fl,]$p
    participant = which(ixes==j)
    s = sign(layerRes[layerRes$layer==layer&layerRes$fl==fl,]$t)
    bg=rgb(s==(-1),0,s==1)
    if (layer==1) {
      points(rep(participant,2),c(-1,100),type='l',col='black')
    }
    points(participant,layer,pch=15,cex=abs(t),col=gray((s+1)/3))
  }
}


