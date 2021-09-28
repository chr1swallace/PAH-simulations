library(magrittr)
library(data.table)

ntests=list(equal.high=rep(10000,7),
            equal.medium=rep(1000,7),
            equal.small=rep(100,7),
            skew.high=c(round(rep(20000/6,6)),50000),
            skew.medium=c(round(rep(2000/6,6)),5000),
            skew.small=c(round(rep(200/6,6)),500),
            data=c(Fig1_B=7,Fig1_Tfh=7,Fig1_Treg=5,Fig2=4,Fig4B=19,Fig4C=20,Fig5=177))
sapply(ntests,sum)

## library(microbenchmark)
## thr=seq(0.1,0.5,length.out=10000)
## f1=function(thr)
##   runif(length(thr)) < thr
## f2=function(thr)
##   sapply(thr, function(a) sample(c(0,1),size=1,prob=c(1-a,a)))
## microbenchmark(f1(thr),
##                f2(thr))

m=as.data.table(ntests[c("skew.small","equal.small","data")])
setnames(m, sub(".small","",names(m)))
m
m=melt(m)
m[,x:=as.character(1:.N),by="variable"]
m[variable=="data",x:=names(ntests$data)]
theme_set(theme_bw())
ggplot(m, aes(x=x,y=value)) +
  geom_col(fill="grey",col="black") +
  facet_wrap(~variable,scales="free_x") +
  theme(legend.position="top",
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  labs(x="stratum",y="stratum size") +
  ggtitle("Number of samples in each stratum in small hypothetical scenarios and observed data",
          sub="Medium and large scenarios have the same shape as small, with 10-fold and 100-fold larger")
ggsave("sfdr-scenarios.png",height=3,width=8)

#' summarise using sens and spec
summ=function(fdr,alt,thr=0.05) {
  c(sens=mean(fdr[alt] < thr),
    spec=1-mean(fdr[!alt] < thr),
    fdr=mean(!alt[fdr < 0.05]))
}
#' simulate
sim=function(altexp=0.1, imbalance=c("random","concentrated"),
             scenario=names(ntests),
             prob.alt=0.3) {
  scenario=match.arg(scenario)
  group=rep(1:length(ntests[[scenario]]),times=ntests[[scenario]])
  imbalance=match.arg(imbalance)
  thr.alt=switch(imbalance,
                 random=prob.alt,
                 concentrated=seq(prob.alt+0.2,prob.alt-0.2,length.out=length(group)))
  alt=runif(length(group)) < thr.alt
  alt.sd= tapply(alt,group,mean) %>% sd
  ## table(group,alt)
  p=exp(-rexp(length(group), ifelse(alt, altexp, 1)))
  ## hist(p[alt])
  ## hist(p[!alt])
  fdr=p.adjust(p,method="BH")
  sfdr=tapply(p, group, p.adjust, method="BH") %>% unlist()
  ## plot(fdr,sfdr)
  data.table(scenario=scenario,imbalance=imbalance, alt.sd, altexp=altexp, prob.alt=prob.alt, ss=c("sens","spec","fdr"), fdr=summ(fdr,alt), sfdr=summ(sfdr,alt))
}

set.seed(42)
results=lapply(names(ntests), function(scen) {
  lapply(rep(c("random","concentrated"),each=1000), function(im)
    sim(imbalance=im,scenario=scen)) %>% rbindlist()
}) %>% rbindlist()
head(results)


library(ggplot2)
m=results[,.(fdr=mean(fdr),sfdr=mean(sfdr),fdr.se=sd(fdr)/sqrt(.N),sfdr.se=sd(sfdr)/sqrt(.N)),by=c("scenario", "imbalance","ss")] %>%
  melt(., c("scenario","imbalance","ss"), measure.vars=list(mean=c("fdr","sfdr"), se=c("fdr.se","sfdr.se"))) #
m[,variable:=c("standard","stratified")[variable]]
m[,scenario:=factor(scenario,levels=names(ntests))]
## just show two classes for simplicity
m[ss=="sens",ss:="Sensitivity"]
m[ss=="fdr",ss:="FDR"]

library(seaborn)
theme_set(theme_bw())
ggplot(m, aes(x=scenario, y=mean,ymin=mean-1.96*se,ymax=mean+1.96*se,col=variable))+
  geom_hline(aes(yintercept=mean+1.96*se),data=m[scenario=="equal.high" & variable=="standard"],linetype="dashed",col="grey") +
  geom_hline(aes(yintercept=mean-1.96*se),data=m[scenario=="equal.high" & variable=="standard"],linetype="dashed",col="grey") +
  geom_pointrange(position = position_dodge(width=.6)) +
  scale_colour_seaborn("FDR approach") +
  labs(y="mean + 95% CI",x="Hypotheses scenario") +
  theme(legend.position="top",
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  facet_grid(ss~imbalance,scales="free")
ggsave("sfdr-results.png",height=8,width=8)

results[,median(alt.sd),by="imbalance"]

rowMeans(results)
