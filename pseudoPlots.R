#!/usr/bin/env Rscript

library(airr)
library(tidyverse)
library(ggplot2)

sj <- function(df) { 
  df %>% 
    group_by(j_call) %>% 
    summarise(freq=length(j_call)/length(df$j_call), repertoire=repertoire[1])
}

sd <- function(df) { 
  df %>% 
    group_by(d_call) %>% 
    summarise(freq=length(d_call)/length(df$d_call), repertoire=repertoire[1])
}


pseudoFreq <- tibble()
repFreq    <- tibble()
jStats     <- tibble()
dStats     <- tibble()
cdrLengths <- tibble()
shmStats   <- tibble()


for( donor in Sys.glob("D*") ) {
  
  setwd(donor)
  for( tissue in Sys.glob("*") ){
    
    for( file in Sys.glob( paste0(tissue, "/FW1/parsedRearrangements.tsv") ) ) {
      
      print(paste(donor,tissue,sep="/"))
      
      r <- read_rearrangement( file )
      r <- r %>% separate( v_call, "[*]", into=c("v_call") )
      r <- r %>% separate( j_call, "[*]", into=c("j_call") )
      r <- r %>% separate( d_call, "[*]", into=c("d_call") )
      
      pvOnly <- r[r$gene_source=="pseudogene",]
      
      pseudoFreq <- rbind( pseudoFreq,
                           pvOnly %>% 
                             group_by(v_call) %>% 
                             summarise(freq=length(v_call)/length(r$v_call), donor=donor, tissue=tissue) )
      
      r$repertoire[ r$gene_source == "pseudogene"] <- "pseudogene"
      repFreq <- rbind( repFreq,
                        r %>% group_by(repertoire) %>%
                          summarise(freq=length(repertoire)/length(r$repertoire), donor=donor, tissue=tissue, source="dna") )
      
      if ( tissue == "Spleen" ) {
        rr <- split(r,r$repertoire)
        
        rj <- lapply(rr, sj)
        for( rept in rj ) { jStats <- rbind(jStats, rept)}
        
        rd <- lapply(rr, sd)
        for( rept in rd ) { dStats <- rbind(dStats, rept)}
        
        cdrLengths <- rbind( cdrLengths, r[,c("junction_length","repertoire")])
        shmStats   <- rbind( shmStats,   r[,c("v_identity", "repertoire")])
        
      }
    }
  }
  
  setwd("..")
  
}

#get RNA data
rnaRep <- tibble()
for( donor in c("40512", "Pt44") ) {
  r <- read_rearrangement( paste0("../2019201_pseudogenes/",donor,"/parsedRearrangements.tsv") )
  r$repertoire[ r$gene_source == "pseudogene"] <- "pseudogene"
  rnaRep <- rbind( rnaRep,
                    r %>% group_by(repertoire) %>%
                      summarise(freq=length(repertoire)/length(r$repertoire), donor=donor, tissue="Spleen", source="H") )
}
for( chain in c("K","L") ) {
  r <- read_rearrangement( paste0("../2019201_pseudogenes/",donor,"/parsedRearrangements.tsv") )
  r$repertoire[ r$gene_source == "pseudogene"] <- "pseudogene"
  rnaRep <- rbind( rnaRep,
                    r %>% group_by(repertoire) %>%
                      summarise(freq=length(repertoire)/length(r$repertoire), donor="40512", tissue="Spleen", source=chain) )
}

#get divergence data
div<-read.table("divergence.txt")
names(div)<-c("v_call","divergence")
usage<-pseudoFreq %>% group_by(v_call) %>% summarise(mean=mean(freq))
mm<-merge(usage,div)
corr<-cor.test(mm$divergence,mm$mean)

###############
#    PLOTS    #
###############

psi_freq <- ggplot( pseudoFreq[pseudoFreq$tissue=="Spleen",] ) + 
  geom_boxplot( aes(x=v_call,y=freq) ) +
  scale_x_discrete(name="VH psuedogene") +
  scale_y_continuous(name="Frequency in VDJ recombinations") +
  theme_bw(base_size=20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),#, vjust = 0.5),
        plot.margin = margin(5.5,5.5,7.5,5.5) )

ggsave("plots/pseudoV_usage.png", psi_freq, w=12)

########

repFreq<-rbind(repFreq,c("pseudogene",NA,"D207","Jejunum"))
repFreq$freq<-as.numeric(repFreq$freq)

tissue_usage <- ggplot( repFreq[repFreq$donor %in% c("D181","D207") & repFreq$repertoire=="pseudogene",] ) + 
  geom_line( aes(x=tissue,y=freq,group=donor,color=donor) ) +
  geom_point( aes(x=tissue,y=freq,color=donor) ) +
  scale_y_continuous(limits = c(0,.01)) +
  labs( x="Tissue", y=expression(paste("Frequency of ",V[p],"DJ recombinations")) ) +
  scale_colour_manual(values = c("red","blue")) +
  theme_bw(base_size=20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),#, vjust = 0.5),
        plot.margin = margin(5.5,5.5,7.5,5.5), legend.position = c(.9,.9), legend.title = element_blank() )

ggsave("plots/tissue_usage.png", tissue_usage, w=10)

#########

expressed <- ggplot( repFreq[repFreq$tissue=="Spleen",] ) + 
  geom_boxplot( aes(x=repertoire,y=freq,group=repertoire,color=donor,shape=source) ) +
  geom_point( data=rnaRep, aes(x=repertoire,y=freq,color=donor,shape=source), size=4, position = "jitter" ) +
  geom_hline( yintercept = 5/7, linetype = 2, color = "#d95f02" ) +
  geom_text( data=data.frame(x=c(2.5),y=c(.8),label=c("theoretical frequncy of functional rearrangements")), aes(x=x,y=y,label=label), color="#d95f02" ) +
  scale_y_continuous(trans="log10") +
  labs( x="Repertoire", y="Frequency" ) +
  scale_colour_manual(values = c("black","#1b9e77","#7570b3"), breaks = c("pooled","40512","Pt44"),
                      limits = c("pooled","40512","Pt44"), labels = c("pooled","40512","Pt44"),
                      guide=guide_legend( override.aes = list(stroke=c(3,3,3), size=c(0,0,0))) ) +
  scale_shape_manual(breaks=c("dna","H","K","L"), values=c(0,15,16,17), 
                     labels=c(expression(gDNA[H]),expression(mRNA[H]),expression(mRNA[K]), expression(mRNA[L])), 
                     guide=guide_legend( override.aes = list(stroke=c(0,5,5,5),size=c(.4,0,0,0))) ) +
  theme_bw(base_size=20)# +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1),#, vjust = 0.5),
#        plot.margin = margin(5.5,5.5,7.5,5.5), legend.position = c(.9,.9), legend.title = element_blank() )

ggsave("plots/expressed.png", expressed, w=10)

########

jgene_usage <- ggplot( jStats ) + 
  geom_boxplot( aes(x=j_call,y=freq,group=interaction(j_call,repertoire),color=repertoire) ) +
  labs( x=expression(paste(J[H]," gene")), y="Frequency" ) +
  scale_colour_manual(values = c("#1b9e77","#d95f02","#7570b3")) +
  theme_bw(base_size=20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),#, vjust = 0.5),
        plot.margin = margin(5.5,5.5,7.5,5.5), legend.position = c(.15,.85), legend.title = element_blank() )

ggsave("plots/jgene_usage.png", jgene_usage, h=6, w=8)

#########

dd<-subset(dStats, !grepl("OR",d_call) & !is.na(d_call) )
dd<- dd %>% separate(d_call,"-",into=c(NA,"pos"),remove=FALSE)
dd$pos <- as.numeric(dd$pos)
dd$d_call <- factor(dd$d_call, levels=unique(dd$d_call[order(dd$pos)]))

dgene_usage <- ggplot( dd ) + 
  geom_boxplot( aes(x=d_call,y=freq,group=interaction(d_call,repertoire),color=repertoire) ) +
  labs( x=expression(paste(D[H]," gene")), y="Frequency" ) +
  scale_colour_manual(values = c("#1b9e77","#d95f02","#7570b3")) +
  theme_bw(base_size=20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),#, vjust = 0.5),
        plot.margin = margin(5.5,5.5,7.5,5.5), legend.position = c(.35,.85), legend.title = element_blank() )

ggsave("plots/dgene_usage.png", dgene_usage, h=6, w=12)

########

junct_length <- ggplot(cdrLengths) +
  geom_density( aes( group=repertoire, x=junction_length, fill=repertoire ), alpha=.4 ) +
  scale_x_continuous( limits = c(0,150), name="junction length (nt)" ) +
  scale_fill_manual(values = c("#1b9e77","#d95f02","#7570b3")) +
  theme_bw(base_size=20) +
  theme( legend.position = c(.85,.85), legend.title = element_blank() )

ggsave("plots/junct_length.png", junct_length, h=6, w=8)

########

shm_density <- ggplot(shmStats) +
  geom_density( aes( group=repertoire, x=1-v_identity, fill=repertoire ), alpha=.4 ) +
  scale_x_continuous( limits = c(0,.4), name=expression(paste("SHM (% ",V[H]," nt)")) ) +
  scale_fill_manual(values = c("#1b9e77","#d95f02","#7570b3")) +
  theme_bw(base_size=20) +
  theme( legend.position = c(.85,.85), legend.title = element_blank() )

ggsave("plots/shm_density.png", shm_density, h=6, w=8)

########

div_vs_use <- ggplot(mm) +
  geom_point( aes( x=divergence, y=mean ) ) +
  geom_text( data=data.frame(x=c(27),y=c(.00087),
		label=c(sprintf("r=%.2f\np=%.2f",corr$estimate[['cor']],corr$p.value))),
 	     aes(x=x,y=y,label=label), size=6 ) +
  labs( x="% nt divergence from nearest functional gene", y="Mean usage frequency") +
  theme_bw(base_size=20) 

ggsave("plots/div_use.png", div_vs_use, h=6, w=8)

