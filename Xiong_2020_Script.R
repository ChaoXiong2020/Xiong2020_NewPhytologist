

### R process ###
rm(list=ls()) 
#########
library(vegan)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(scales)
library(grid)
library(dplyr)
library(tidyr)
library(splines)
library(nlme)
library(lmerTest)
library(lme4)
library(pbkrtest)
library(ade4)
library(lsmeans)

####### Alpha diversity #######

#### Import data ####
setwd("F:/Desktop/r/chaoxiong/crop_microbiome/alpha/")
alpha=read.csv(file="alpha.csv",header=T,check.names=FALSE,sep=",")
head(alpha)
#### Plot Figure  ###
p=ggplot(alpha, aes(x=niche, y=shannon, fill=niche,col=niche))+ 
  geom_boxplot(width=0.5,position=position_dodge(0.6),alpha=0.8)+geom_jitter(width = 0.2,alpha=0.8,size=1)
p=p+scale_colour_manual(values =c("gray33","#999999","#E69F00","#009E73","#D55E00","#56B4E9"))
p=p+scale_fill_manual(values =c("gray33","#999999","#E69F00","#009E73","#D55E00","#56B4E9"))
p=p+theme(axis.title.x = element_text(face="bold", colour="black", size=6),
            axis.title.y = element_text(face="bold", colour="black", size=6),
            axis.text.x = element_text(face="bold", colour="black", size=6),
            axis.text.y  = element_text(face="bold", colour="black", size=6))
p
#####

####### Beta diversity #######

#### Import data ####
setwd("F:/Desktop/r/chaoxiong/crop_microbiome/beta/")
beta=read.csv(file="alpha.csv",header=T,check.names=FALSE,sep=",")
head(beta)
#### Plot Figure  ###
p=ggplot(beta, aes(x=xvar, y=yvar, color=niche,shape=crop)) + geom_point(size=1.2,alpha=0.7)
p=p+scale_shape_manual(name="Crop",breaks = c('maize','wheat','barley'),values = c(16,17,15))
p=p+scale_colour_manual(name="niche",values =c("#009E73","#56B4E9","#E69F00","#D55E00", "#999999","gray33"))
p=p+theme(axis.title.x = element_text(face="bold", colour="black", size=6),
            axis.title.y = element_text(face="bold", colour="black", size=6),
            axis.text.x = element_text(face="bold", colour="black", size=6),
            axis.text.y  = element_text(face="bold", colour="black", size=6))
p
#####

####### Community composition #######

#### Import data ####
setwd("F:/Desktop/r/chaoxiong/crop_microbiome/composition/")
sp=read.csv(file="composition.csv",header=T,check.names=FALSE,sep=",")
head(sp)
df <- sp %>% gather("item",value,-1) %>% 
bind_cols(data.frame(item_id=rep(1:12,each=684))) 
head(df)
df$item = factor(df$item, levels=c('Gammaproteobacteria','Alphaproteobacteria','Betaproteobacteria',
                                   'Deltaproteobacteria','Actinobacteria','Bacteroidetes',
                                   'Firmicutes','Saccharibacteria','Gemmatimonadetes',
                                   'Chloroflexi','Acidobacteria','Other'))
#### Plot Figure  ###								   
p=ggplot(df,aes(sample,value))+geom_bar(aes(fill=item),stat = "identity",position="fill",width=0.8)
p=p+scale_fill_manual(name=NULL,values = c("dark green","green","springgreen3","cadetblue3","purple",
                                           "orange","violetred","brown","red",'#D55E00', "blue4","grey50","tan","salmon",
                                           "dark blue","yellow","black"))
p=p+theme(axis.title.x = element_text(face="bold", colour="black", size=6),
            axis.title.y = element_text(face="bold", colour="black", size=6),
            axis.text.x = element_text(face="bold", colour="black", size=6),
            axis.text.y  = element_text(face="bold", colour="black", size=6))
p=p + scale_y_continuous(labels=percent)
p
#####

####### The linear mixed model #######

#### Import data ####
setwd("F:/Desktop/r/chaoxiong/crop_microbiome/LMM/")
sp=read.csv(file="LMM.csv",header=T,check.names=FALSE,sep=",")
head(sp)
### Analysis ####
div.fit1=lmer(shannon~Host niche + Crop species + Site + Fertilization practice + (1|block), REML=TRUE,data=sp)
Anova(div.fit1,type=2,ddf="lme4",test="F")
r.squaredGLMM(div.fit1)
#####

####### PERMANOVA #######

#### Import data ####
setwd("F:/Desktop/r/chaoxiong/crop_microbiome/PERMANOVA/")
sp.dist=read.table(file="weighted_unifrac_zotu_table_norm_maize.txt",header=T,check.names=FALSE,sep="")
map=read.table(file="map_maize.txt",header=T,row.names=1,check.names=FALSE,sep="")
head(sp.dist)
head(map)
### Analysis ####
adonis_maize=adonis(sp.dist~niche+site+Treatment,data =map,permutations=999)

### Qiime process ###

### CSS normalize process ###
normalize_table.py -i ../qiime_output/zotu_table.biom -a CSS -o ../qiime_output/zotu_table_norm.biom

### Alpha_diversity ###
alpha_diversity.py -i zotu_table_even.biom -m chao1,shannon,simpson,simpson_e,heip_e,goods_coverage,observed_species,PD_whole_tree -t zotus.tree -o alpha.txt

### Beta_diversity ###
beta_diversity.py -i ../qiime_output/zotu_table_norm.biom -m weighted_unifrac -t ../output/zotus.tree -o ../qiime_output/beta_diversity

### Beta_diversity_nmds ###
nmds.py -i ../qiime_output/beta_diversity -o ../qiime_output/nmds

### Summarize_Community composition ###
summarize_taxa.py -i ../qiime_output/zotu_table.biom -o ../qiime_output/tax_summary_a -L 1,2,3,4,5,6,7 -a 
summarize_taxa.py -i ../qiime_output/zotu_table.biom -o ../qiime_output/tax_summary_r -L 1,2,3,4,5,6,7

### Source-tracking analysis based on SourceTracker (v1.0) ###
R --slave --vanilla --args -i otus.txt -m metadata.txt -o sourcetracker_out < sourcetracker_for_qiime.r

### Function prediction based on FAPROTAX ###
python collapse_table.py -i zotu_table.txt -o fun_table.txt -g FAPROTAX.txt -c '#' -d 'taxonomy' --omit_columns 0 --column_names_are_in last_comment_line -r report.txt -n columns_after_collapsing -v


