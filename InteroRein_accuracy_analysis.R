library(psycho)
#https://neuropsychology.github.io/psycho.R/2018/03/29/SDT.html

##########################################################################################
##########################################################################################

##########################################################################################
#RECOGONITION
#get percentage of signal detection aggregate by number of correctl/wrong trials
#recognition accuracy by congruency and cycle
rec.acc.df_all <- aggregate(list(rec.acc.df$trial_count,rec.acc.df$trial_total_count),by=(list(rec.acc.df$congruency, rec.acc.df$cardiac_cycle, rec.acc.df$image_category, rec.acc.df$sigDetection, rec.acc.df$subject, rec.acc.df$subject_condition)), FUN=sum)
colnames(rec.acc.df_all) <- c("congruency", "cardiac_cycle", "image_category",  "sigDetection", "subject", "subject_condition", "trial_count", "trial_total_count")
row.names(rec.acc.df_all) <- NULL
#set up two groups of subject
rec.acc.df_subset <- rec.acc.df_all[rec.acc.df_all$subject_condition=='diastole_face+systole_scene',]
rec.acc.df_subset <- rec.acc.df_all[rec.acc.df_all$subject_condition=='diastole_scene+systole_face',]
rec.acc.df_subset <- rec.acc.df_all
# # to deal with 100 hit rate and 0 false alarm, we're using the option to: "add 0.5 to both the number of hits and the number of false alarms, and add 1 to both the number of signal trials and the number of noise trials; dubbed the loglinear approach (Hautus, 1995)"
# # https://stats.stackexchange.com/questions/134779/d-prime-with-100-hit-rate-probability-and-0-false-alarm-probability
# for (n in 1:nrow(rec.acc.df_subset)){
#   if (rec.acc.df_subset$trial_count[n]==0 & (n%%4==1 | n%%4==3)){
#     rec.acc.df_subset$trial_count[n] <- rec.acc.df_subset$trial_count[n]+0.5
#     rec.acc.df_subset$trial_count[n+1] <- rec.acc.df_subset$trial_count[n+1]+0.5
#     rec.acc.df_subset$trial_total_count[n] <- rec.acc.df_subset$trial_total_count[n]+1
#     rec.acc.df_subset$trial_total_count[n+1] <- rec.acc.df_subset$trial_total_count[n+1]+1
#   }else if (rec.acc.df_subset$trial_count[n]==0 & (n%%4==2 | n%%4==0)){
#     rec.acc.df_subset$trial_count[n] <- rec.acc.df_subset$trial_count[n]+0.5
#     rec.acc.df_subset$trial_count[n-1] <- rec.acc.df_subset$trial_count[n-1]+0.5
#     rec.acc.df_subset$trial_total_count[n] <- rec.acc.df_subset$trial_total_count[n]+1
#     rec.acc.df_subset$trial_total_count[n-1] <- rec.acc.df_subset$trial_total_count[n-1]+1
#   }
#   
#   if (rec.acc.df_subset$trial_count[n]==rec.acc.df_subset$trial_total_count[n] & (n%%4==1 | n%%4==3)){
#     rec.acc.df_subset$trial_count[n] <- rec.acc.df_subset$trial_count[n]+0.5
#     rec.acc.df_subset$trial_count[n+1] <- rec.acc.df_subset$trial_count[n+1]+0.5
#     rec.acc.df_subset$trial_total_count[n] <- rec.acc.df_subset$trial_total_count[n]+1
#     rec.acc.df_subset$trial_total_count[n+1] <- rec.acc.df_subset$trial_total_count[n+1]+1
#   }else if (rec.acc.df_subset$trial_count[n]==rec.acc.df_subset$trial_total_count[n] & (n%%4==2 | n%%4==0)){
#     rec.acc.df_subset$trial_count[n] <- rec.acc.df_subset$trial_count[n]+0.5
#     rec.acc.df_subset$trial_count[n-1] <- rec.acc.df_subset$trial_count[n-1]+0.5
#     rec.acc.df_subset$trial_total_count[n] <- rec.acc.df_subset$trial_total_count[n]+1
#     rec.acc.df_subset$trial_total_count[n-1] <- rec.acc.df_subset$trial_total_count[n-1]+1
#   }
# }
rec.acc.df_subset$percentage <- rec.acc.df_subset$trial_count/rec.acc.df_subset$trial_total_count

################################
# cardiac cycle and congrency  
################################
#plot d' or A'
dprim_df <-  data.frame()
for (a in unique(rec.acc.df_subset$congruency)){
  for (b in unique(rec.acc.df_subset$cardiac_cycle)){
    dprime <- c()
    aprime <- c()
    for (x in seq_along(unique(rec.acc.df_subset$subject))){
      df_temp <- rec.acc.df_subset[rec.acc.df_subset$subject==unique(rec.acc.df_subset$subject)[x],]
      dprime[x] <- dprime(n_hit = sum(df_temp$trial_count[df_temp$sigDetection=='hit' & df_temp$congruency==a & df_temp$cardiac_cycle==b]), 
               n_fa = sum(df_temp$trial_count[df_temp$sigDetection=='false alarm' & df_temp$congruency==a & df_temp$cardiac_cycle==b]), 
               n_miss = sum(df_temp$trial_count[df_temp$sigDetection=='miss' & df_temp$congruency==a & df_temp$cardiac_cycle==b]), 
               n_cr = sum(df_temp$trial_count[df_temp$sigDetection=='correct rejection' & df_temp$congruency==a & df_temp$cardiac_cycle==b]))$dprime
      aprime[x] <- dprime(n_hit = sum(df_temp$trial_count[df_temp$sigDetection=='hit' & df_temp$congruency==a & df_temp$cardiac_cycle==b]), 
                          n_fa = sum(df_temp$trial_count[df_temp$sigDetection=='false alarm' & df_temp$congruency==a & df_temp$cardiac_cycle==b]), 
                          n_miss = sum(df_temp$trial_count[df_temp$sigDetection=='miss' & df_temp$congruency==a & df_temp$cardiac_cycle==b]), 
                          n_cr = sum(df_temp$trial_count[df_temp$sigDetection=='correct rejection' & df_temp$congruency==a & df_temp$cardiac_cycle==b]))$aprime
      }
    dprim_tem <- data.frame("d.prime"=dprime,"a.prime"=aprime,"cardiac_cycle"=rep(b,length(dprime)),"congruency"=rep(a,length(dprime)), "subject"=unique(rec.acc.df_subset$subject))
    dprim_df <- rbind(dprim_df,dprim_tem)
  }
}
ggplot(dprim_df, aes(x=interaction(cardiac_cycle, congruency),y=d.prime)) + 
  stat_summary(fun=mean,position=position_dodge(width=0.90), alpha = 0.3, geom="bar", size = 1.5)+
  stat_summary(fun.data=mean_se, geom="errorbar",width=0.1,  size = 0.5)+
  geom_violin(trim=TRUE, alpha = 0.05, size=0.02)+
  geom_line(data = dprim_df, aes(x = interaction(cardiac_cycle, congruency), y = d.prime, group=subject, color=subject), alpha = 0.3, size = 0.5)+
  geom_point(data = dprim_df, aes(x = interaction(cardiac_cycle, congruency), y = d.prime, group=subject, color=subject), alpha = 0.3, size = 1)+
  # facet_wrap(congruency~cardiac_cycle)+
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=16,face="bold"),
        strip.text = element_text(size = 13)) 
summary(aov(d.prime ~ cardiac_cycle * congruency + Error(subject/(cardiac_cycle * congruency)), data=dprim_df))
summary(aov(a.prime ~ cardiac_cycle * congruency + Error(subject/(cardiac_cycle * congruency)), data=dprim_df), na.action=na.omit)
# t.test(dprim_df[dprim_df$congruency=='congruent' & dprim_df$cardiac_cycle=='diastole',]$d.prime,dprim_df[dprim_df$congruency=='congruent' & dprim_df$cardiac_cycle=='systole',]$d.prime, paired=TRUE)
# t.test(dprim_df[dprim_df$congruency=='incongruent' & dprim_df$cardiac_cycle=='diastole',]$d.prime,dprim_df[dprim_df$congruency=='incongruent' & dprim_df$cardiac_cycle=='systole',]$d.prime, paired=TRUE)
# dprime(n_hit = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='hit' & rec.acc.df$image_category=='scene']), n_fa = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='false alarm' & rec.acc.df$image_category=='scene']), n_miss = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='miss' & rec.acc.df$image_category=='scene']), n_cr = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='correct rejection' & rec.acc.df$image_category=='scene']))
# dprime(n_hit = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='hit' & rec.acc.df$image_category=='face']), n_fa = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='false alarm' & rec.acc.df$image_category=='face']), n_miss = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='miss' & rec.acc.df$image_category=='face']), n_cr = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='correct rejection' & rec.acc.df$image_category=='face']))

#plot false alarm and hit rate
dprime_labels <- aggregate(list(dprim_df$a.prime, dprim_df$d.prime),by=(list(dprim_df$cardiac_cycle, dprim_df$congruency)), FUN=mean, na.rm=TRUE)
colnames(dprime_labels) <- c("cardiac_cycle", "congruency", "a'", "d'")
dprime_labels$sigDetection<-'hit'
dprime_labels$`a'`<-paste0("A:", round(dprime_labels$`a'`, digits = 3))
dprime_labels$`d'`<-paste0("d:", round(dprime_labels$`d'`, digits = 3))
ggplot(rec.acc.df_subset[rec.acc.df_subset$sigDetection=='hit' | rec.acc.df_subset$sigDetection=='false alarm',], aes(x=sigDetection,y=percentage, fill=sigDetection)) + 
  stat_summary(fun=mean,position=position_dodge(width=0.90), alpha = 0.3, geom="bar", size = 1.5)+
  stat_summary(fun.data=mean_se, geom="errorbar",width=0.1,  size = 0.5)+
  geom_violin(trim=TRUE, alpha = 0.05, size=0.02)+
  geom_line(data = rec.acc.df_subset[rec.acc.df_subset$sigDetection=='hit' | rec.acc.df_subset$sigDetection=='false alarm',], aes(x = sigDetection, y = percentage, group=subject, color=subject), alpha = 0.3, size = 0.5)+
  geom_point(data = rec.acc.df_subset[rec.acc.df_subset$sigDetection=='hit' | rec.acc.df_subset$sigDetection=='false alarm',], aes(x = sigDetection, y = percentage, group=subject, color=subject), alpha = 0.3, size = 1)+
  facet_grid(congruency~cardiac_cycle)+
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=16,face="bold"),
        strip.text = element_text(size = 13)) +
  geom_text(aes(sigDetection, 0.9, label = `a'`), data = dprime_labels, parse = TRUE, hjust = 0)+
  geom_text(aes(sigDetection, 1, label = `d'`), data = dprime_labels, parse = TRUE, hjust = 0)

################################
# cardiac cycle and image category
################################
#plot d' or A'
dprim_df <-  data.frame()
for (a in unique(rec.acc.df_subset$image_category)){
  for (b in unique(rec.acc.df_subset$cardiac_cycle)){
    dprime <- c()
    aprime <- c()
    for (x in seq_along(unique(rec.acc.df_subset$subject))){
      df_temp <- rec.acc.df_subset[rec.acc.df_subset$subject==unique(rec.acc.df_subset$subject)[x],]
      dprime[x] <- dprime(n_hit = sum(df_temp$trial_count[df_temp$sigDetection=='hit' & df_temp$image_category==a & df_temp$cardiac_cycle==b]), 
                          n_fa = sum(df_temp$trial_count[df_temp$sigDetection=='false alarm' & df_temp$image_category==a & df_temp$cardiac_cycle==b]), 
                          n_miss = sum(df_temp$trial_count[df_temp$sigDetection=='miss' & df_temp$image_category==a & df_temp$cardiac_cycle==b]), 
                          n_cr = sum(df_temp$trial_count[df_temp$sigDetection=='correct rejection' & df_temp$image_category==a & df_temp$cardiac_cycle==b]))$dprime
      aprime[x] <- dprime(n_hit = sum(df_temp$trial_count[df_temp$sigDetection=='hit' & df_temp$image_category==a & df_temp$cardiac_cycle==b]), 
                          n_fa = sum(df_temp$trial_count[df_temp$sigDetection=='false alarm' & df_temp$image_category==a & df_temp$cardiac_cycle==b]), 
                          n_miss = sum(df_temp$trial_count[df_temp$sigDetection=='miss' & df_temp$image_category==a & df_temp$cardiac_cycle==b]), 
                          n_cr = sum(df_temp$trial_count[df_temp$sigDetection=='correct rejection' & df_temp$image_category==a & df_temp$cardiac_cycle==b]))$aprime
    }
    dprim_tem <- data.frame("d.prime"=dprime,"a.prime"=aprime,"cardiac_cycle"=rep(b,length(dprime)),"image_category"=rep(a,length(dprime)), "subject"=unique(rec.acc.df_subset$subject))
    dprim_df <- rbind(dprim_df,dprim_tem)
  }
}
ggplot(dprim_df, aes(x=interaction(cardiac_cycle, image_category),y=d.prime)) + 
  stat_summary(fun=mean,position=position_dodge(width=0.90), alpha = 0.3, geom="bar", size = 1.5)+
  stat_summary(fun.data=mean_se, geom="errorbar",width=0.1,  size = 0.5)+
  geom_violin(trim=TRUE, alpha = 0.05, size=0.02)+
  geom_line(data = dprim_df, aes(x = interaction(cardiac_cycle, image_category), y = d.prime, group=subject, color=subject), alpha = 0.3, size = 0.5)+
  geom_point(data = dprim_df, aes(x = interaction(cardiac_cycle, image_category), y = d.prime, group=subject, color=subject), alpha = 0.3, size = 1)+
  # facet_wrap(image_category~cardiac_cycle)+
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=16,face="bold"),
        strip.text = element_text(size = 13)) 
summary(aov(d.prime ~ cardiac_cycle * image_category + Error(subject/(cardiac_cycle * image_category)), data=dprim_df))
summary(aov(a.prime ~ cardiac_cycle * image_category + Error(subject/(cardiac_cycle * image_category)), data=dprim_df), na.action=na.omit)
# t.test(dprim_df[dprim_df$image_category=='congruent' & dprim_df$cardiac_cycle=='diastole',]$d.prime,dprim_df[dprim_df$image_category=='congruent' & dprim_df$cardiac_cycle=='systole',]$d.prime, paired=TRUE)
# t.test(dprim_df[dprim_df$image_category=='incongruent' & dprim_df$cardiac_cycle=='diastole',]$d.prime,dprim_df[dprim_df$image_category=='incongruent' & dprim_df$cardiac_cycle=='systole',]$d.prime, paired=TRUE)
# dprime(n_hit = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='hit' & rec.acc.df$image_category=='scene']), n_fa = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='false alarm' & rec.acc.df$image_category=='scene']), n_miss = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='miss' & rec.acc.df$image_category=='scene']), n_cr = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='correct rejection' & rec.acc.df$image_category=='scene']))
# dprime(n_hit = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='hit' & rec.acc.df$image_category=='face']), n_fa = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='false alarm' & rec.acc.df$image_category=='face']), n_miss = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='miss' & rec.acc.df$image_category=='face']), n_cr = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='correct rejection' & rec.acc.df$image_category=='face']))

#plot false alarm and hit rate
dprime_labels <- aggregate(list(dprim_df$a.prime, dprim_df$d.prime),by=(list(dprim_df$cardiac_cycle, dprim_df$image_category)), FUN=mean, na.rm=TRUE)
colnames(dprime_labels) <- c("cardiac_cycle", "image_category", "a'", "d'")
dprime_labels$sigDetection<-'hit'
dprime_labels$`a'`<-paste0("A:", round(dprime_labels$`a'`, digits = 3))
dprime_labels$`d'`<-paste0("d:", round(dprime_labels$`d'`, digits = 3))
ggplot(rec.acc.df_subset[rec.acc.df_subset$sigDetection=='hit' | rec.acc.df_subset$sigDetection=='false alarm',], aes(x=sigDetection,y=percentage, fill=sigDetection)) + 
  stat_summary(fun=mean,position=position_dodge(width=0.90), alpha = 0.3, geom="bar", size = 1.5)+
  stat_summary(fun.data=mean_se, geom="errorbar",width=0.1,  size = 0.5)+
  geom_violin(trim=TRUE, alpha = 0.05, size=0.02)+
  geom_line(data = rec.acc.df_subset[rec.acc.df_subset$sigDetection=='hit' | rec.acc.df_subset$sigDetection=='false alarm',], aes(x = sigDetection, y = percentage, group=subject, color=subject), alpha = 0.3, size = 0.5)+
  geom_point(data = rec.acc.df_subset[rec.acc.df_subset$sigDetection=='hit' | rec.acc.df_subset$sigDetection=='false alarm',], aes(x = sigDetection, y = percentage, group=subject, color=subject), alpha = 0.3, size = 1)+
  facet_grid(image_category~cardiac_cycle)+
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=16,face="bold"),
        strip.text = element_text(size = 13)) +
  geom_text(aes(sigDetection, 0.9, label = `a'`), data = dprime_labels, parse = TRUE, hjust = 0)+
  geom_text(aes(sigDetection, 1, label = `d'`), data = dprime_labels, parse = TRUE, hjust = 0)

################################
# cardiac cycle and image category and congruency
################################
#plot d' or A'
dprim_df <-  data.frame()
for (a in unique(rec.acc.df_subset$image_category)){
  for (b in unique(rec.acc.df_subset$cardiac_cycle)){
    for (c in unique(rec.acc.df_subset$congruency)){
      dprime <- c()
      aprime <- c()
      subject_condition_temp <- c()
      for (x in seq_along(unique(rec.acc.df_subset$subject))){
        df_temp <- rec.acc.df_subset[rec.acc.df_subset$subject==unique(rec.acc.df_subset$subject)[x],]
        dprime[x] <- dprime(n_hit = sum(df_temp$trial_count[df_temp$sigDetection=='hit' & df_temp$image_category==a & df_temp$cardiac_cycle==b & df_temp$congruency==c]), 
                            n_fa = sum(df_temp$trial_count[df_temp$sigDetection=='false alarm' & df_temp$image_category==a & df_temp$cardiac_cycle==b & df_temp$congruency==c]), 
                            n_miss = sum(df_temp$trial_count[df_temp$sigDetection=='miss' & df_temp$image_category==a & df_temp$cardiac_cycle==b & df_temp$congruency==c]), 
                            n_cr = sum(df_temp$trial_count[df_temp$sigDetection=='correct rejection' & df_temp$image_category==a & df_temp$cardiac_cycle==b & df_temp$congruency==c]))$dprime
        aprime[x] <- dprime(n_hit = sum(df_temp$trial_count[df_temp$sigDetection=='hit' & df_temp$image_category==a & df_temp$cardiac_cycle==b & df_temp$congruency==c]), 
                            n_fa = sum(df_temp$trial_count[df_temp$sigDetection=='false alarm' & df_temp$image_category==a & df_temp$cardiac_cycle==b & df_temp$congruency==c]), 
                            n_miss = sum(df_temp$trial_count[df_temp$sigDetection=='miss' & df_temp$image_category==a & df_temp$cardiac_cycle==b & df_temp$congruency==c]), 
                            n_cr = sum(df_temp$trial_count[df_temp$sigDetection=='correct rejection' & df_temp$image_category==a & df_temp$cardiac_cycle==b & df_temp$congruency==c]))$aprime
        subject_condition_temp[x] <- rec.acc.df_subset$subject_condition[rec.acc.df_subset$subject==unique(rec.acc.df_subset$subject)[x]][1]
        }
      dprim_tem <- data.frame("d.prime"=dprime,"a.prime"=aprime,
                              "cardiac_cycle"=rep(b,length(dprime)),"image_category"=rep(a,length(dprime)),"congruency"=rep(c,length(dprime)), 
                              "subject"=unique(rec.acc.df_subset$subject),"subject_condition"=subject_condition_temp)
      dprim_df <- rbind(dprim_df,dprim_tem)
    }
  }
}
dprim_df <- na.omit(dprim_df)
ggplot(dprim_df, aes(x=interaction(cardiac_cycle, image_category),y=d.prime)) + 
  stat_summary(fun=mean,position=position_dodge(width=0.90), alpha = 0.3, geom="bar", size = 1.5)+
  stat_summary(fun.data=mean_se, geom="errorbar",width=0.1,  size = 0.5)+
  geom_violin(trim=TRUE, alpha = 0.05, size=0.02)+
  geom_line(data = dprim_df, aes(x = interaction(cardiac_cycle, image_category), y = d.prime, group=subject, color=subject), alpha = 0.3, size = 0.5)+
  geom_point(data = dprim_df, aes(x = interaction(cardiac_cycle, image_category), y = d.prime, group=subject, color=subject), alpha = 0.3, size = 1)+
  facet_wrap(subject_condition~congruency,drop = TRUE)+
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=16,face="bold"),
        strip.text = element_text(size = 13)) 
summary(aov(d.prime ~ cardiac_cycle * image_category + Error(subject/(cardiac_cycle * image_category)), data=dprim_df))
summary(aov(a.prime ~ cardiac_cycle * image_category + Error(subject/(cardiac_cycle * image_category)), data=dprim_df), na.action=na.omit)
# t.test(dprim_df[dprim_df$image_category=='congruent' & dprim_df$cardiac_cycle=='diastole',]$d.prime,dprim_df[dprim_df$image_category=='congruent' & dprim_df$cardiac_cycle=='systole',]$d.prime, paired=TRUE)
# t.test(dprim_df[dprim_df$image_category=='incongruent' & dprim_df$cardiac_cycle=='diastole',]$d.prime,dprim_df[dprim_df$image_category=='incongruent' & dprim_df$cardiac_cycle=='systole',]$d.prime, paired=TRUE)
# dprime(n_hit = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='hit' & rec.acc.df$image_category=='scene']), n_fa = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='false alarm' & rec.acc.df$image_category=='scene']), n_miss = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='miss' & rec.acc.df$image_category=='scene']), n_cr = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='correct rejection' & rec.acc.df$image_category=='scene']))
# dprime(n_hit = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='hit' & rec.acc.df$image_category=='face']), n_fa = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='false alarm' & rec.acc.df$image_category=='face']), n_miss = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='miss' & rec.acc.df$image_category=='face']), n_cr = sum(rec.acc.df$trial_count[rec.acc.df$sigDetection=='correct rejection' & rec.acc.df$image_category=='face']))

#plot false alarm and hit rate
dprime_labels <- aggregate(list(dprim_df$a.prime, dprim_df$d.prime),by=(list(dprim_df$cardiac_cycle, dprim_df$image_category)), FUN=mean, na.rm=TRUE)
colnames(dprime_labels) <- c("cardiac_cycle", "image_category", "a'", "d'")
dprime_labels$sigDetection<-'hit'
dprime_labels$`a'`<-paste0("A:", round(dprime_labels$`a'`, digits = 3))
dprime_labels$`d'`<-paste0("d:", round(dprime_labels$`d'`, digits = 3))
ggplot(rec.acc.df_subset[rec.acc.df_subset$sigDetection=='hit' | rec.acc.df_subset$sigDetection=='false alarm',], aes(x=sigDetection,y=percentage, fill=sigDetection)) + 
  stat_summary(fun=mean,position=position_dodge(width=0.90), alpha = 0.3, geom="bar", size = 1.5)+
  stat_summary(fun.data=mean_se, geom="errorbar",width=0.1,  size = 0.5)+
  geom_violin(trim=TRUE, alpha = 0.05, size=0.02)+
  geom_line(data = rec.acc.df_subset[rec.acc.df_subset$sigDetection=='hit' | rec.acc.df_subset$sigDetection=='false alarm',], aes(x = sigDetection, y = percentage, group=subject, color=subject), alpha = 0.3, size = 0.5)+
  geom_point(data = rec.acc.df_subset[rec.acc.df_subset$sigDetection=='hit' | rec.acc.df_subset$sigDetection=='false alarm',], aes(x = sigDetection, y = percentage, group=subject, color=subject), alpha = 0.3, size = 1)+
  facet_grid(image_category~cardiac_cycle)+
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=16,face="bold"),
        strip.text = element_text(size = 13)) +
  geom_text(aes(sigDetection, 0.9, label = `a'`), data = dprime_labels, parse = TRUE, hjust = 0)+
  geom_text(aes(sigDetection, 1, label = `d'`), data = dprime_labels, parse = TRUE, hjust = 0)

