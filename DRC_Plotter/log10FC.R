#Log10 FC Pipeline

#Select COnstruct for Normalization
constructFC=sub_df[4, 'construct'] #FOr testing
constructFC_ec50=sub_df[sub_df$construct==constructFC, 'ec50']

#Subset Dataframe
FC_df <- sub_df %>% select(ligand, construct, n, 'log10FC_ec50'=ec50, 'log10FC_ec50sd'=ec50_sd)


#log10 Transform Dataframe
FC_df$log10FC_ec50<-(log10(FC_df$log10FC_ec50)-log10(constructFC_ec50))/log10(constructFC_ec50)
FC_df$log10FC_ec50sd<-log10(FC_df$log10FC_ec50sd)
FC_sd_cols=c('FC50_upper', 'FC50_lower')
FC_df[, FC_sd_cols]<-NA
FC_df$FC50_upper<-FC_df$log10FC_ec50+FC_df$log10FC_ec50sd
FC_df$FC50_lower<-FC_df$log10FC_ec50-FC_df$log10FC_ec50sd

#Calculate Bounds
FC_max<-abs(max(ceiling(FC_df$FC50_upper)))
FC_min<-abs(min(floor(FC_df$FC50_lower)))
FC_bounds=max(FC_max, FC_min)

#Filter-out Ligand that data is Normalized to:
FC_norm <- FC_df %>% filter(construct!=constructFC)

#Build Barplot
limits <- aes(ymax = FC50_upper,
              ymin = FC50_lower)

ggplot(data = FC_norm,
       aes(x = construct, y = log10FC_ec50)
) +
  geom_bar(
    stat = "identity",
    aes(fill = log10FC_ec50),
    position = position_dodge(width = 0.9)
  ) +
  ylab("Fold Change (log_10)") +
  xlab("Construct") +
  geom_errorbar(
    limits,
    width = 0.15,
    position = position_dodge(width = 0.9)
  ) +
  coord_flip(ylim=c(-FC_bounds, FC_bounds)) +
  ggtitle(paste0("Log10 Fold Change of EC50 \n (Normalized to ", constructFC, ")"))+
  theme(plot.title = element_text(hjust = 0.5))
