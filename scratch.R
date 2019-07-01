#lseq
library(emdbook)
library(ggplot2)
# library(readxl)
library(scales)
library(RColorBrewer)
library(utils)
library(Hmisc)
library(tidyverse)
library(dplyr)


#HARCODED PATH FOR TESTING
# data_path="C:\\Users\\Vasco Morais\\Documents\\computation\\hill_plotter\\test_data.csv"
#HARCODED PATH FOR TESTING
# data_path="C:\\Users\\Vasco Morais\\Documents\\computation\\hill_plotter\\duplicate_test.csv"
#SUM STAT PATH
data_path='Z:\\Shared\\Company\\Confidential\\Research\\Team IonFlux\\Dose Response Curves\\summary_stats_master_df.csv'

#Import Data
p_df <- read.csv(data_path, header = T, stringsAsFactors=FALSE, strip.white = TRUE)

#FUNCTIONS
#Load and preclean-up data
import_data = function(p_df) {
  
  #Uses defaults if summar_stats_master_df
  if ('Dose_Sum' %in% colnames(p_df)) {
    
        #Import specific columns
    p_df <- p_df %>%
      select(Exp, Plate_Type, Ligand, construct, Key, bottom, top, ec50, hill, success_rate)
    #Fix column Headings
    colnames(p_df) <- tolower(colnames(p_df))
    col_names <- colnames(p_df)
    
    #Variable for Automatic Sum_stat DF
    df_type='auto'
    
  } else {
    
    #FIX COL HEADING
    col_names<-c('key', 'construct',	'bottom',	'top',	'ec50',	'hill')
    colnames(p_df) <- col_names
    
    #Variable for Automatic Sum_stat DF
    df_type='manual'
    
  }
  
  #remove rows without data
  p_df=p_df[col_names][complete.cases(p_df[col_names]),]
  
  #Remove duplicate rows
  p_df<-unique(p_df)
  
  # Removes rows where 'top' bounds is under 1000
  p_df=p_df %>% filter(top>1000)
  
  #Convert list of ligands to lowercase
  p_df$key<-sapply(p_df$key, tolower)
  
  #get rid of problematic naming elements
  p_df<-mutate_if(p_df, 
                  is.character, 
                  str_replace_all, pattern = "-", replacement = "_")
  #ONLY IF A7_4
  #change 'a7_4' to 'pcoda102'
  if ('a7_4' %in% p_df$construct) {
    p_df[p_df['construct']=='a7_4',]['construct']='CODA102'
  }
  
  return(list(p_df, df_type))
}

#Geometric Mean Function
"geometric.mean" <- 
  function(x,na.rm=TRUE)
  { 
    exp(mean(log(x),na.rm=na.rm))
  }

#GEOMETRIC STDEV
geom_sd=function(gmean, n, A) {
  exp(sqrt(sum(log(A/gmean)**2)/n))
}


#Assign Ligand based off of lig-key
ligand_ind <- function(key) {
  if (startsWith(key, 'ach')) {
    lig='Acetylcholine'
  } else if (startsWith(key, 'fac')) {
    lig='Facinicline'
  } else if (startsWith(key, 'nic')) {
    lig='Nicotine'
  } else if (startsWith(key, 'gly')) {
    lig='Glycine'
  } else if (startsWith(key, 'ivm')) {
    lig='Ivermectin'
  } else if (startsWith(key, 'azd')) {
    lig='AZD-0328'
  } else if (startsWith(key, 'abt')) {
    lig='ABT-126'
  } else if (startsWith(key, 'tc')) {
    lig='TC-6987'
  } else {
    lig='other'
  }
  return(lig)
}


#Iterates through dataframe renaming any duplicate constructs that are present within a specific ligand
rename_duplicate_constructs <- function(df) {
  #Identifies duplicate constructs for a given ligand type
  duplicate_key <- as.data.frame(df %>% group_by(ligand) %>% count(construct) %>% filter(nn>1) %>% select(ligand, construct))
  
  #Append lig_key name to duplicate constructs
  if (nrow(duplicate_key)>0) {
    for (r in 1:nrow(duplicate_key)) {
      l=duplicate_key[r, 'ligand'] #ligand
      c=duplicate_key[r, 'construct'] #construct
      
      keys=df[(df['ligand']==l & df['construct']==c),]$key
      
      #Changes name of construct to following format: "Construct [Key]"
      for (key in keys) {
        df[(df['ligand']==l & df['construct']==c & df['key']==key),]$construct <- paste0(c, ' [', toupper(gsub('_', '.', key)), ']')
      }
    }
  }
  return(df)
}


#Hill Equation
logistic_array = function(x, beta) {
  return(beta[1] + (beta[2]-beta[1])/(1 + (x/beta[3])^(-beta[4])))
}


# Logistic function for dose-response curve, compact coefficients
# Parameters
# -----------
# beta = list/array of length 4
# beta[1] = bottom
# beta[2] = top
# beta[3] = ec50
# beta[4] = hill_coefficient
# x = Ligand concentration
# Returns
# -----------
# f(x): predicted ligand current




#Calculate Bounds & Error Bars
roundUp=function(x) {10^ceiling(log10(x))}
roundDown=function(x) {10^floor(log10(x))}

bounds=function(y, beta) {
  
  inverse_eq=function(y, beta) {beta[3]*exp(((beta[2]-beta[1])/(-beta[4]*(y*(beta[2])-beta[1])))+(1/beta[4]))}
  bounds=inverse_eq(y, beta)
  return(bounds)
  
}

plot_bounds=function(df) {
  bounds_lower=array()
  bounds_upper=array()
  for(j in 1:nrow(df)) {
    
    beta=as.numeric(df[j,c('bottom', 'top', 'ec50', 'hill')])
    
    bounds_lower=c(bounds_lower, bounds(0.25, beta))
    bounds_upper=c(bounds_upper, roundUp(bounds(0.9, beta))*10)
    
  }
  
  bounds_lower=bounds_lower[!is.na(bounds_lower)]
  bounds_lower=roundDown(min(bounds_lower))
  
  bounds_upper=bounds_upper[!is.na(bounds_upper)]
  bounds_upper=roundUp(max(bounds_upper))
  
  bounds=c(bounds_lower, bounds_upper)
  
  return(bounds)
}


#Append population (n) to names
construct_n <- function(j, df) {
  n=df[j, 'n'] #population
  c=df[j, 'construct'] #construct
  
  construct_pop=paste0(c, ' (n=', n, ')')
  
  return(construct_pop)
}






##### PIPELINE BEGINS HERE #####


imported_data=import_data(p_df)

# Cleaned-Up Dataframe
p_df=imported_data[[1]]

# Dataframe type (auto or manual)
df_type=imported_data[[2]]

#Calculate Population Counts
n_df<-as.data.frame(
  p_df %>%
    group_by(key) %>%
    count(construct)
  )

#Merge population counts onto DF
df<-inner_join(p_df, n_df, by=c('key', 'construct'))

#Calculate Geometric Mean
if (df_type=='auto') {
  df <- as.data.frame(
    df %>%
      group_by(exp, plate_type, ligand, key, construct, n) %>%
      summarise(bottom=geometric.mean(bottom), top=geometric.mean(top), ec50=geometric.mean(ec50), hill=geometric.mean(hill), success_rate=mean(success_rate))
  )
} else {
  df <- as.data.frame(
    df %>%
      group_by(key, construct, n) %>%
      summarise(bottom=geometric.mean(bottom), top=geometric.mean(top), ec50=geometric.mean(ec50), hill=geometric.mean(hill))
  ) 
}



#Creates empty columns to be filled with stats
if (df_type=='auto') {
    ec50_Cols=c('ec50_sd', 'ec50_upper', 'ec50_lower', 'ec50_y_Norm')
    df[, ec50_Cols]<-NA
} else {
    ec50_Cols_manual=c('ec50_sd', 'ec50_upper', 'ec50_lower', 'ec50_y_Norm', 'ligand')
    df[, ec50_Cols_manual]<-NA
}



#Calculate Geometric Stats and Assign Ligand name
for (r in 1:nrow(df)) {
  k=df[r, 'key'] #key
  c=df[r, 'construct'] #construct
  
  #Gets list of all ec50 values for lig-receptor pair
  A=p_df[p_df['key']==k & p_df['construct']==c,]$ec50 #ec50 value
  
  #n Population value
  n=df[r, 'n']
  
  #Mean ec50
  m=df[r, 'ec50']
  
  #Calculate geometric Standard Deviation and save into df
  df[r, 'ec50_sd']<-geom_sd(m, n, A)
  
  #Calculate Geometric Standard Deviation Upper and Lower Bounds
  if (df[r, 'n']>1) {
    #Calculates Upper Bounds
    df[r, 'ec50_upper'] <- df[r, 'ec50']*df[r, 'ec50_sd']
    #Calculates Lower Bounds
    df[r, 'ec50_lower'] <- df[r, 'ec50']/df[r, 'ec50_sd']
  }
  
  #Calculate ec50_y_Norm
  beta=as.numeric(df[r, c('bottom', 'top', 'ec50', 'hill')])
  
  #IF NOT NORMALIZED (REACTIVE)
  # df[r, 'ec50_y_Norm']<-logistic_array(m, beta)
  
  #IF NORMALIZED (REACTIVE)
  df[r, 'ec50_y_Norm']<-logistic_array(m, beta)/beta[2]*100
  
  if (df_type=='manual') {
    #Assign Ligand name based off key
    df[r, 'ligand']<-ligand_ind(df[r, 'key'])
  }
}


df=rename_duplicate_constructs(df)


#Assign colors
if (length(unique(df$construct)) <= 8) {
  construct.colors=brewer.pal(length(unique(df$construct)), "Dark2") #Dark2 works for 8 constructs in total, Set3 is one color larger
  names(construct.colors)<-levels(as.factor(unique(df$construct)))#Loads in constructs and assigns to color
}


# #Applies population (n) to end of construct name
# for (i in 1:nrow(df)) {
#   con_n=construct_n(i, df)
#   df[i, ]$construct<-con_n
# }



#Calculate Hill Equation (REACTIVE - Normalization)
for (lig in unique(df$ligand)) {
  print(lig)
  ligsub_df<-df %>% group_by(ligand) %>% filter(ligand==lig)
  ligsub_df<-as.data.frame(ligsub_df)

  hill_eq=data.frame()

  for(j in 1:nrow(ligsub_df)) {

    beta=as.numeric(ligsub_df[j,c('bottom', 'top', 'ec50', 'hill')])

    x_range=lseq(plot_bounds(ligsub_df)[1], plot_bounds(ligsub_df)[2], 500)

    drug_curve=logistic_array(x = x_range, beta)/beta[2]*100

    #Normalized to pcoda71
    # p71_top=eval(parse(text=i))[eval(parse(text=i))['construct']=='pcoda71', 'top']
    # drug_curve=logistic_array(x = x_range, beta)/p71_top

    construct=rep(ligsub_df[j,'construct'], 500)

    hill_eq=rbind(hill_eq, as.data.frame(cbind.data.frame(x_range, drug_curve, construct)))

    colnames(hill_eq)=c('x_range', 'drug_curve', 'construct')
  }


  if (length(unique(df$construct)) <= 8) {
  drug_plot=ggplot(hill_eq, aes(x=x_range, y=drug_curve)) +
    geom_line(size=2, alpha=0.5, aes(group=construct, color=construct)) +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = 'b') +
    scale_color_manual(values=construct.colors)
  } else {
    drug_plot=ggplot(hill_eq, aes(x=x_range, y=drug_curve)) +
      geom_line(size=2, alpha=0.5, aes(group=construct, color=construct)) +
      scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      annotation_logticks(sides = 'b') #+
      # scale_color_manual(values=construct.colors)
  }
  
  #Produce DRC Plot
  drug_plot_aes=drug_plot +
    xlab(bquote(.(lig) ~ "["*mu*M*"]")) +
    ylab('Normalized Peak Current \n (% of Imax)') +
    theme_grey(base_size = 22) +
    theme(legend.title = element_blank()) +
    ggtitle(paste(ligsub_df[j,2])) + #overwritten below
    ggtitle(label = NULL)

  # Subset rows that only have a standard dev greater than 1
  errorbar_df=ligsub_df[!is.na(ligsub_df['ec50_upper'] | ligsub_df['ec50_lower']),]

  drug_plot_stat= drug_plot_aes +
    geom_errorbarh(data=errorbar_df, mapping=aes(y=ec50_y_Norm, x=ec50, xmin=ec50_lower, xmax=ec50_upper, color=construct), height=7, size=1.5, alpha=0.5) +
    geom_point(data=ligsub_df, mapping=aes(group=construct, y=ec50_y_Norm, x=ec50, color=construct), size=5, shape=20, alpha=0.3, na.rm = TRUE)

  print(drug_plot_stat)

}

