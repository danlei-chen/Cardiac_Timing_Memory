library(reshape2)
library(ggplot2)
library(effsize)
library(tidyr)
library(psycho)
library(lme4)
library(nlme)

subj_list <- c("002","003","004","005","012","014","015","017","018","019","020","021","022","023","024","026","027","028","029","030","031")
data_dir <- '/Volumes/IASL/Studies/InteroceptiveReinstatement/subjData/'

enc.df <- data.frame()
enc.acc_all <- c()
for (subj in subj_list){
  print(subj)
  
  #encoding
  subj.enc.df <- read.csv(file = paste0(data_dir,'subj-',subj,'_InteroRein_encoding_output/subj-',subj,'_Encoding_Response_by_Trial.csv'))
  subj.enc.df$block_category<- 'face'
  subj.enc.df$block_category[(subj.enc.df$trial_label=='indoor' | subj.enc.df$trial_label=='outdoor')] <- 'scene'
  subj.enc.df$subject <- subj
  if (!("missing.heartbeat" %in% colnames(subj.enc.df))){
    subj.enc.df$missing.heartbeat <- ""
  }
  enc.df <- rbind(enc.df, subj.enc.df)
  
  enc.acc <- aggregate(subj.enc.df$trial_accuracy,by=(list(subj.enc.df$trial_accuracy)), FUN=length)
  enc.acc_all <- c(enc.acc_all, enc.acc$x[enc.acc$Group.1==1]/(enc.acc$x[enc.acc$Group.1==1]+enc.acc$x[enc.acc$Group.1==0]))
  print(paste0('enc: ',enc.acc_all[length(enc.acc_all)]))
}
print(paste0('grand mean of all subject:', mean(enc.acc_all)))

###################
#build encoding df that contains variables that we want
enc.df$trial_accuracy[enc.df$trial_accuracy=='NaN']<-0
#create a new array ind filled with trial_accuracy
#use dcast to use value.var column to aggregation function, length, in this case. And returns counts/length value to the column after ~
enc.acc.df.all <- dcast(transform(enc.df, ind='trial_accuracy'), trial_accuracy + block_num + subject + cardiac_cycle + block_category~ind,
                        value.var='trial_accuracy', length, drop=FALSE)
colnames(enc.acc.df.all) <- c("trial_accuracy","block","subject","cardiac_cycle","block_category","trial_accuracy_count")

#set the accuracy df as half of the orginal and make accuracy 
enc.acc.df <- enc.acc.df.all[enc.acc.df.all$trial_accuracy==1,]
colnames(enc.acc.df)[6] <- 'correct_trial_count'
enc.acc.df$total_trial_count <- enc.acc.df$correct_trial_count + enc.acc.df.all$trial_accuracy_count[enc.acc.df.all$trial_accuracy==0]
#make accuracy using num correct divided by total count
enc.acc.df$accuracy <- enc.acc.df$correct_trial_count/enc.acc.df$total_trial_count

enc.acc.df <- enc.acc.df[!is.na(enc.acc.df$accuracy),]
enc.acc.df <- subset(enc.acc.df, select = -c(trial_accuracy) )
enc.acc.df$categoryAndCycle <- paste0(enc.acc.df$cardiac_cycle, '_', enc.acc.df$block_category)
enc.acc.df$subject_condition <- ifelse(enc.acc.df$categoryAndCycle == 'diastole_scene' | enc.acc.df$categoryAndCycle == 'systole_face', 'diastole_scene+systole_face',
                                      ifelse(enc.acc.df$categoryAndCycle == 'diastole_face' | enc.acc.df$categoryAndCycle == 'systole_scene', 'diastole_face+systole_scene','others'))

################################################
#subject/block sanity check
#encoding accuracy by block
ggplot(enc.acc.df,aes(block,accuracy))+
  geom_point(data = enc.acc.df, aes(block,accuracy, group=subject, color=subject),alpha = 0.3, size = 1.5)+
  stat_summary(fun.data=mean_se, geom="ribbon", alpha = 0.15, colour = NA)+
  geom_line(data = enc.acc.df, aes(block,accuracy, group=subject, color=subject), alpha = 0.3, size = 0.5)+
  stat_summary(fun.data=mean_se, geom="line", size = 0.5)+
  stat_summary(fun=mean, geom="point", size = 1.5)+ 
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=16,face="bold")) + 
  ggtitle("Encoding accuracy by block")

################################################
#Group condition comparison, hopefully not different
#aggregate cross subject for plottings
enc.acc.df.temp <- aggregate(enc.acc.df$accuracy, by=(list(enc.acc.df$subject, enc.acc.df$subject_condition)), FUN=mean, na.rm=TRUE)
colnames(enc.acc.df.temp) <- c("subject","subject_condition","accuracy")
ggplot(enc.acc.df.temp,aes(subject_condition,accuracy), group=subject)+
  stat_summary(fun=mean, geom="bar", size = 1.5, alpha=0.2)+ 
  geom_jitter(data = enc.acc.df.temp, aes(subject_condition,accuracy, group=subject, color=subject),alpha = 0.5, size = 1.5, width = 0.05, height = 0)+
  stat_summary(fun.data=mean_se, geom="errorbar",width=0.1,  size = 0.5)+
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=16,face="bold")) + 
  ggtitle("Group condition") +
  coord_cartesian(ylim=c(min(enc.acc.df.temp$accuracy), NA))
t.test(enc.acc.df.temp$accuracy[enc.acc.df.temp$subject_condition=='diastole_scene+systole_face'], enc.acc.df.temp$accuracy[enc.acc.df.temp$subject_condition=='diastole_face+systole_scene'], paired=FALSE)

################################################
#full comparison on each condition
#with block
ggplot(enc.acc.df,aes(block,accuracy, group=interaction(block_category,cardiac_cycle), col=interaction(block_category,cardiac_cycle), fill = interaction(block_category,cardiac_cycle)))+
  stat_summary(fun.data=mean_se, geom="ribbon", alpha = 0.15, colour = NA)+
  stat_summary(fun.data=mean_se, geom="line", size = 0.5)+
  # geom_point(data = enc.acc.df, aes(block,accuracy, group=subject, color=subject),alpha = 0.3, size = 1.5)+
  # geom_line(data = enc.acc.df, aes(block,accuracy, group=subject, color=subject), alpha = 0.3, size = 0.5)+
  stat_summary(fun=mean, geom="point", size = 1.5)+ 
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=15),
        plot.title = element_text(size=18,face="bold"),
        strip.text = element_text(size=15)) + 
  ggtitle("Block category by block") + 
  facet_wrap(~subject_condition)
#CANNOT use repeated measure on cardiac_cycle https://www.r-bloggers.com/how-to-do-repeated-measures-anovas-in-r/
#repeated measure ANOVA on each group 
summary(aov(accuracy ~ categoryAndCycle * block + Error(subject/(categoryAndCycle * block)), data=enc.acc.df[enc.acc.df$subject_condition=='diastole_face+systole_scene',]))
summary(aov(accuracy ~ categoryAndCycle * block + Error(subject/(categoryAndCycle * block)), data=enc.acc.df[enc.acc.df$subject_condition=='diastole_scene+systole_face',]))
#mixed effect ANOVA using subject_condition as random effect
summary(aov(accuracy ~ subject_condition * categoryAndCycle * block + Error(subject/(categoryAndCycle * block)), data=enc.acc.df))
summary(aov(accuracy ~ subject_condition * block_category * cardiac_cycle * block + Error(subject/(block_category * cardiac_cycle * block)), data=enc.acc.df))
# anova(fit2<-lmer(accuracy ~ cardiac_cycle * block_category * block + (1|subject), data=enc.acc.df,na.action = na.exclude))

################################################
#collapse across block
enc.acc.df_noblock <- aggregate(enc.acc.df$accuracy, by=(list(enc.acc.df$subject, enc.acc.df$block_category, enc.acc.df$cardiac_cycle, enc.acc.df$categoryAndCycle, enc.acc.df$subject_condition)), FUN=mean, na.rm=TRUE)
colnames(enc.acc.df_noblock) <- c("subject", "block_category", "cardiac_cycle", "categoryAndCycle", "subject_condition", "accuracy")
enc.acc.df_noblock.mean <- aggregate(enc.acc.df$accuracy, by=(list(enc.acc.df$subject_condition)), FUN=mean, na.rm=TRUE)
colnames(enc.acc.df_noblock.mean) <- c("subject_condition", "mean_accuracy")
                                  
ggplot(enc.acc.df_noblock,aes(categoryAndCycle,accuracy),group=categoryAndCycle, color=subject)+ 
  geom_violin(trim=TRUE, alpha = 0.1, size=0.1)+
  geom_point(data = enc.acc.df_noblock, aes(categoryAndCycle,accuracy, group=subject, color=subject), alpha = 0.25, size = 1)+
  geom_line(data = enc.acc.df_noblock, aes(categoryAndCycle,accuracy, group=subject, color=subject), alpha = 0.25, size = 0.5)+
  stat_summary(fun=mean,position=position_dodge(width=0.90), alpha = 0.3, geom="bar", size = 1.5)+
  stat_summary(fun.data=mean_se,position=position_dodge(width=0.90), geom="errorbar",width=0.1,  size = 0.5)+ 
  geom_hline(data= enc.acc.df_noblock.mean, aes(yintercept=mean_accuracy), linetype='dotted') + 
  # stat_smooth(method="rq", formula=accuracy~1, se=FALSE)+
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=15),
        plot.title = element_text(size=18,face="bold"),
        strip.text = element_text(size=15)) + 
  ggtitle("Encoding accuracy by cardiac cycle") + coord_cartesian(ylim=c(min(enc.acc.df_noblock$accuracy), NA)) + 
  facet_wrap(~subject_condition, scales="free_x")
t.test(enc.acc.df_noblock$accuracy[enc.acc.df_noblock$subject_condition=='diastole_scene+systole_face' & enc.acc.df_noblock$categoryAndCycle=='diastole_scene'], enc.acc.df_noblock$accuracy[enc.acc.df_noblock$subject_condition=='diastole_scene+systole_face' & enc.acc.df_noblock$categoryAndCycle=='systole_face'], paired=FALSE)
t.test(enc.acc.df_noblock$accuracy[enc.acc.df_noblock$subject_condition=='diastole_face+systole_scene' & enc.acc.df_noblock$categoryAndCycle=='diastole_face'], enc.acc.df_noblock$accuracy[enc.acc.df_noblock$subject_condition=='diastole_face+systole_scene' & enc.acc.df_noblock$categoryAndCycle=='systole_scene'], paired=FALSE)



