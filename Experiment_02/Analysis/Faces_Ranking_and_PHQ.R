rm(list=ls())

library(lm.beta)
library(ggplot2)
library(rstudioapi)

# Fix terrible plotting on Windows
trace(grDevices:::png, quote({if (missing(type) && missing(antialias)) {
  type <- "cairo-png"
  antialias <- "subpixel"}}), print = FALSE)

# Get current path
script_path = rstudioapi::getActiveDocumentContext()$path
pathSplit=strsplit(script_path, "/")
pathSplit=pathSplit[[1]]
main_path=paste0(pathSplit[1:(length(pathSplit)-3)],"/",collapse="")


## Sample II - Replication
#For iMac
path=paste0(main_path,"Boost_faces_emotional_II/Output/")
subjects=c(101:102,104:114,116:118,120:122,124:125,127:139,142:152,155:157,159:166,168:181)
# exclude:
# 103 - participates in a previous CAT experiment (and had high false alarm in this experiment) 
# 123, 126, 140, 154, 167 - Training: false alarm
# 158 - Training: misses
# 115, 158 - Training: minimal ladder
# 119, 141, 153 - quit in the middle

## Sample I
path=paste0(main_path,"Boost_faces_emotional/Output/")
subjects=c(113,301,303:304,306:307,309:312,314:316,318:326,328:331,333:346,348:351,353:355,357:359) 
#subjects = c(113,132,301:312,314:326,328:331,333:346,348:359)
# exclude:
# 132, 302, 305 - Training: Minimal Ladder
# 308, 317 - Training: false alarm
# 356 - Binary ranking
# (technical)
# % 313 - was coded as 113
# % 327 - code crashed at BR
# % 332 - was coded as 132 
# % 347 - quit in the middle of the experiment

filelist=c()
for (s in subjects){
  # These files were created with the MATAB code: 'CreateItemRanking.m'
  filelist=c(filelist,Sys.glob(paste(path, "*",s,"*ItemRankingResults_combined*.txt",sep="")))
}
Data=c()
for (s in subjects){
  tmp_file=Sys.glob(paste(path, "*",s,"*ItemRankingResults_combined*.txt",sep=""))
  tmp_data=read.table(tmp_file,header=T,na.strings=c(999,999000))
  differences = tmp_data$Rank[1:40] - tmp_data$Rank[41:80]
  Happy_bias = mean(differences)/sd(differences)
  SubjectID = as.character(tmp_data$Subject[1])
  Data=rbind(Data,data.frame(SubjectID,Happy_bias))
}

PHQ_data=read.table(Sys.glob(paste(path, "/PHQ_Scores.txt",sep="")),header=T)
Data = merge(Data,PHQ_data,by='SubjectID')

Data$PHQ_ordinal=rank(Data$PHQ)
Data$Happy_bias_ordinal=rank(Data$Happy_bias)

model_linear=lm(PHQ ~ Happy_bias,data=Data)
summary(lm.beta(model_linear))

dev.new()
ggplot(data=Data,aes(Happy_bias,PHQ,color=PHQ)) +
  geom_point(shape = 16, size = 5, alpha= 0.8) +
  geom_vline(xintercept=0,linetype=2) +
  geom_smooth(method='lm') +
  xlab("Initial Preferences - Bias towards Happy Faces") +
  theme_bw()

model_ordinal=lm(PHQ_ordinal ~ Happy_bias_ordinal,data=Data)
summary(lm.beta(model_ordinal))

dev.new()
ggplot(data=Data,aes(Happy_bias_ordinal,PHQ_ordinal,color=PHQ)) +
  geom_point(shape = 16, size = 5, alpha= 0.8) +
  geom_smooth(method='lm') +
  xlab("Initial Preferences - Bias towards Happy Faces (Rank)") +
  ylab("PHQ (Rank)") +
  theme_bw()

# Plot PHQ distribution - 2 groups ----

Data$PHQ_binary_n = 1
Data$PHQ_binary_n[Data$PHQ>=9] = 2
Data$PHQ_binary = factor(Data$PHQ_binary_n, labels = c("Low", "High"))

ggplot(Data,aes(PHQ_binary,PHQ,color=PHQ_binary)) +
  set.seed(0) + geom_point(shape = 16, size = 4, alpha= 0.6,position = position_jitter(w = 0.5, h = 0)) +
  set.seed(0) + geom_point(shape = 21, size = 4, alpha= 1,position = position_jitter(w = 0.5, h = 0), color = "black") +
  geom_hline(yintercept = 9, linetype=2) +
  geom_hline(yintercept = 5, linetype=2) +
  xlab("Experimental group") + ylab("PHQ")+
  theme(text  = element_text(size = 20)) + 
  scale_y_continuous(lim = c(0,22)) +
  theme_bw()

# Plot PHQ distribution - 3 groups ----

Data$PHQ_binary_n = 2
Data$PHQ_binary_n[Data$PHQ<=4] = 1
Data$PHQ_binary_n[Data$PHQ>=10] = 3
Data$PHQ_binary = factor(Data$PHQ_binary_n, labels = c("Low", "Med", "High"))

ggplot(Data,aes(PHQ_binary,PHQ,color=PHQ_binary)) +
  set.seed(0) + geom_point(shape = 16, size = 4, alpha= 0.6,position = position_jitter(w = 0.5, h = 0)) +
  set.seed(0) + geom_point(shape = 21, size = 4, alpha= 1,position = position_jitter(w = 0.5, h = 0), color = "black") +
  geom_hline(yintercept = 9, linetype=2) +
  geom_hline(yintercept = 5, linetype=2) +
  xlab("Experimental group") + ylab("PHQ")+
  theme(text  = element_text(size = 20)) + 
  scale_y_continuous(lim = c(0,22)) +
  theme_bw()

ggplot(Data,aes(PHQ_binary,PHQ,color=PHQ_binary)) +
  geom_point(shape = 21, size = 4, alpha= 1,position = position_jitterdodge(jitter.width = 3)) +
  geom_hline(yintercept = 9, linetype=2) +
  geom_hline(yintercept = 5, linetype=2) +
  xlab("Experimental group") +
  scale_y_continuous(lim = c(0,22)) +
  theme_bw()

ggplot(Data,aes(PHQ_binary,PHQ,color=PHQ_binary)) +
  geom_point(shape = 21, size = 2, alpha= 1,position = position_dodge2(width = 1, padding = 0.5, preserve = "total")) +
  geom_hline(yintercept = 9, linetype=2) +
  geom_hline(yintercept = 5, linetype=2) +
  xlab("Experimental group") +
  scale_y_continuous(lim = c(0,22)) +
  theme_bw()

ggplot(Data,aes(x=PHQ,fill=PHQ_binary)) +
  geom_histogram(binwidth = 1) +
  geom_vline(xintercept = 9, linetype=2) +
  geom_vline(xintercept = 5, linetype=2) +
  xlab("PHQ") +
  scale_y_continuous() +
  theme_bw()

