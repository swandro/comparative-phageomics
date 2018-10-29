##Visualization of hit matrices
library(reshape2)

#Import data
nt.files <- list.files("hit_matrix/nt")

#Format data
nt.df <- NULL
for (file in nt.files){
  openfile <- read.delim(paste("hit_matrix/nt/",file,sep = ''), header = T, sep = '\t',check.names = F)
  name = unlist(strsplit(file, "_hit"))[1]
  openfile[name] <- 100
  melted <- melt(openfile, id.vars = "gene")
  melted$source <- name
  nt.df <- rbind(nt.df, melted)
}

#Fix V12 name
nt.df$variable <- as.character(nt.df$variable)
nt.df$variable[which(nt.df$variable=="V12_RAST")] <- "V12"
nt.df$source[which(nt.df$source=="V12_RAST")] <- "V12"

#Order phages on y axis
levels(factor(nt.df$variable))
nt.df$variable = factor(nt.df$variable, levels = c("Bop","V12","Ben","Bill","CCS1","NC_029009",
                                                   "KY549443","TEX","NC_029026","NC_027335",
                                                   "AP018714","AP018715","Car","Carl","Bob",
                                                   "phief24c_oriented"),ordered=T)


#Plot
for (phage in levels(factor(nt.df$source))){
  P <- ggplot(data = subset(nt.df, source == phage)) +
    geom_tile( aes(x= gene, y=variable, fill = value)) +
    scale_fill_gradient(high="#31a354", low="#e5f5e0",na.value = "white")  +
    scale_y_discrete(expand = c(.005,0)) +  
    labs(x=paste(phage,"genes"), y='', fill="ANI %") +
    theme(text=element_text(size=8),
          axis.text.x = element_blank(),
          axis.text.y= element_text(color="black", size=8),
          panel.background = element_rect(color="black",fil=NA),
          axis.line=element_line(),
          panel.grid.major.x =  element_blank(),
          panel.grid.minor.y= element_line(color="black",size=2),
          legend.key.size=unit(.15,"in"))
  ggsave(filename = paste("plots/genome_conservation/",phage,".png",sep = ''),plot = P,device = "png",height = 2,width=6)
}
       