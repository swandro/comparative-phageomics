##Visualization of hit matrices
library(reshape2)

#Import data
nt.files <- list.files("hit_matrix/nt")

nt.df <- NULL
for (file in nt.files){
  openfile <- read.delim(paste("hit_matrix/nt/",file,sep = ''), header = T, sep = '\t',check.names = F)
  name = unlist(strsplit(file, "_hit"))[1]
  openfile[name] <- 100
  melted <- melt(openfile, id.vars = "gene")
  melted$source <- name
  nt.df <- rbind(nt.df, melted)
}

nt.df$variable <- as.character(nt.df$variable)
nt.df$variable[which(nt.df$variable=="V12_RAST")] <- "V12"
nt.df$source[which(nt.df$source=="V12_RAST")] <- "V12"

levels(factor(nt.df$variable))
nt.df$variable = factor(nt.df$variable, levels = c("Bop","V12","Ben","Bill","CCS1","NC_029009",
                                                   "KY549443","TEX","NC_029026","NC_027335",
                                                   "AP018714","AP018715","Car","Carl","Bob",
                                                   "phief24c_oriented"),ordered=T)
str(nt.df)


ggplot(data = subset(nt.df, source == "NC_029026")) +
  geom_tile( aes(x= gene, y=variable, fill = value)) +
  scale_fill_gradient(high="#31a354", low="#e5f5e0",na.value = "white")  +
  theme(axis.text.x = element_blank())

temp <- subset(nt.df, source == "Car")
subset(temp, is.na(variable)
       
       