#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

#For Pipeline
library(emdbook)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(utils)
library(Hmisc)
library(tidyverse)
library(dplyr)
library(lubridate)
library(readxl)

#For WebApp
library(shiny)
library(shinythemes)
library(DT)
# library(shinyWidgets)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  ###REACTIVE VALUES
  RV_df <- reactiveValues(
    final_df=data.frame(),
    ligands='',
    
    ligkey='',
    constructs='',
    
    exp_range=''
  )
  
  RV_plot <- reactiveValues(
    drc=ggplot(), #DRC Plot
    pop_df=data.frame()
  )
  
  RV_FC <- reactiveValues(
    fc_constructs="",
    fc_plot=ggplot()
  )
  

  #START OF CODE
  observeEvent(c(input$unprocesssed_sumstat_df, input$ligand, input$key, input$construct, input$normalize, input$exp_range, input$run, input$plate_type, input$success_rate, input$geombars, input$constructFC), {
    
    ####FUNCTIONS
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
        
        #Fix exp_range column
        #Removes characters from exp_id
        p_df$exp <- gsub("[^0-9.]", "", p_df$exp)
        #Next two lines change exp_id format form yymmdd to yyyymmdd format (only works for 2018 and 2019 dates that are written in yymmdd format)
        p_df$exp <- gsub("^18", "2018", p_df$exp)
        p_df$exp <- gsub("^19", "2019", p_df$exp)
        #Converts exp_id into yyyy-mm-dd string
        p_df$exp <- ymd(p_df$exp)
        
        #Variable for Automatic Sum_stat DF
        df_type='auto'
        
      } else if ('Timings' %in% colnames(p_df)) {
        
        p_df <- p_df[1:13]
        
        #ColNames
        manColnames<-c('exp', 'construct', 'timings', 'ligand', 'key', 'bottom', 'top', 'ec50', 'hill', 'peakcurrent_75', 'peakcurrent_max', 'fit', 'success_rate')

        #Import and Rename Columns from Manual Summary Table
        colnames(p_df) <- manColnames
        # colnames(p_df)[1] <- 'exp'
        
        # p_df <- select(p_df, exp='exp', timings='Timings', ligand="Ligand", key='Ligand.Key', construct="CODA.Construct",
               # bottom='Bottom', top='Top', ec50='EC50', hill='Hill.Slope',
               # peakcurrent_75="Efficacy..Peak.Current....75.", peakcurrent_max="Efficacy...Peak.Current..max.",
               # fit="Goodness.of.Curve.Fit..Subjective.", success_rate="Success.Rate..good..stable..responsive.sweeps.")
        
        
        #Remove success_rate column (It's not percentage and can't be used to subset)
        p_df <- select(p_df, -success_rate)
        
        
        #Capture column names
        col_names<- colnames(p_df)
        
        #Fix exp_range column
        #Removes characters from exp_id
        p_df$exp <- gsub("[^0-9.]", "", p_df$exp)
        #Next two lines change exp_id format form yymmdd to yyyymmdd format (only works for 2018 and 2019 dates that are written in yymmdd format)
        p_df$exp <- gsub("^18", "2018", p_df$exp)
        p_df$exp <- gsub("^19", "2019", p_df$exp)
        #Converts exp_id into yyyy-mm-dd string
        p_df$exp <- ymd(p_df$exp)
        
        #Convert columns to correct datatypes
        cols.num <- c('bottom', 'top', 'ec50', 'hill', 'peakcurrent_75', 'peakcurrent_max')
        p_df[cols.num] <- sapply(p_df[cols.num], as.numeric)
        cols.int <- 'fit'
        p_df[cols.int] <- sapply(p_df[cols.int], as.integer)
        
        #Variable for Automatic Sum_stat DF
        df_type='manual'
        
      } else {
        
        #FIX COL HEADING
        col_names<-c('key', 'construct',	'bottom',	'top',	'ec50',	'hill')
        colnames(p_df) <- col_names
        
        #Variable for Automatic Sum_stat DF
        df_type='other'
        
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
      duplicate_key <- as.data.frame(df %>% group_by(ligand) %>% count(construct) %>% filter(n>1) %>% select(ligand, construct))
      
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
    
    
    ####START OF PIPELINE
    
    #IMPORTS FILE
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    inFile <- input$unprocesssed_sumstat_df
    
    if (is.null(inFile))
      return(NULL)

    if (endsWith(inFile$name, '.xlsx')) {
      p_df <- read_xlsx(inFile$datapath, sheet = 'Summary Table', progress = TRUE)
      p_df <- as.data.frame(apply(p_df,2,trimws), stringsAsFactors = FALSE)

    } else if (endsWith(inFile$name, '.csv')) {
      p_df<-read.csv(inFile$datapath, header = TRUE, stringsAsFactors=FALSE, strip.white = TRUE)
    } else {
      print('Filetype Incompatible')
    }
      
    
    #Load and preclean-up data
    imported_data=import_data(p_df)
    
    # Cleaned-Up Dataframe
    p_df=imported_data[[1]]
    
    # Dataframe type (auto or manual)
    df_type=imported_data[[2]]
    
    
    #Print initial datframe out for subsetting by exp_range
    if (df_type!='other') {
      RV_df$final_df <- p_df
    }
    
    
    # Exp Range
    # RV_df$exp_range <- unique(sub_df$exp)
    # sub_df<-filter(sub_df, exp>=input$exp_range[1] & exp<=input$exp_range[2])
    
    # Exp Range
    # RV_df$exp_range <- unique(p_df$exp)
    if (is.na(input$exp_range[1])==FALSE & is.na(input$exp_range[2])==FALSE & df_type!='other') {
      #Subsets dataframe by the exp_range chosen
      p_df<-filter(p_df, exp>=input$exp_range[1] & exp<=input$exp_range[2])
      RV_df$final_df <- p_df
    }
    
    
    
    
    
    
    
    
    #PRESUBSET DF BY EXP_RANGE (REQUIRED TO RUN)
    if (df_type=='other' | (is.na(input$exp_range[1])==FALSE & is.na(input$exp_range[2])==FALSE & input$run==TRUE & nrow(p_df)>=1)) {
        
      
        #Calculate Population Counts
        n_df<-as.data.frame(
          p_df %>%
            group_by(key) %>%
            count(construct)
        )
        
        #Merge population counts onto DF
        df<-inner_join(p_df, n_df, by=c('key', 'construct'))
        
        #Calculate Geometric Mean (Calculate either geom-mean or arithmetric mean)
        if (df_type=='auto') {
          df <- as.data.frame(
            df %>%
              group_by(plate_type, ligand, key, construct, n) %>%
              summarise(bottom=geometric.mean(bottom), top=geometric.mean(top), ec50=geometric.mean(ec50), hill=geometric.mean(hill), success_rate=mean(success_rate))
          )
        } else if (df_type=='manual') {
          
          df <- as.data.frame(
          df %>%
            group_by(ligand, key, construct, n) %>%
            summarise(bottom=geometric.mean(bottom), top=geometric.mean(top), ec50=geometric.mean(ec50), hill=geometric.mean(hill))
          )
        } else {
          df <- as.data.frame(
            df %>%
              group_by(key, construct, n) %>%
              summarise(bottom=geometric.mean(bottom), top=geometric.mean(top), ec50=geometric.mean(ec50), hill=geometric.mean(hill))
          ) 
        }
        
        
        #Creates empty columns to be filled with stats
        if (df_type=='auto' | df_type=='manual') {
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
          
          #EC50_y Normalization (REACTIVE)
          if (input$normalize) {
            df[r, 'ec50_y_Norm']<-logistic_array(m, beta)/beta[2]*100
          } else {
            df[r, 'ec50_y_Norm']<-logistic_array(m, beta)
          }
          
          
          if (df_type=='other') {
            #Assign Ligand name based off key
            df[r, 'ligand']<-ligand_ind(df[r, 'key'])
          }
        }
        
        
        df=rename_duplicate_constructs(df)
        
        
        # #Assign colors
        # if (length(unique(df$construct)) <= 8) {
        #   construct.colors=brewer.pal(length(unique(df$construct)), "Dark2") #Dark2 works for 8 constructs in total, Set3 is one color larger
        #   names(construct.colors)<-levels(as.factor(unique(df$construct)))#Loads in constructs and assigns to color
        # }
        
        #Reorganize DF for readibility    
        if (df_type=='manual') {
          df <- df
        } else if (df_type=='auto') {
          df <- df %>%
            select(plate_type, ligand, key, construct, n, bottom, top, ec50, hill, success_rate, ec50_sd, ec50_upper, ec50_lower, ec50_y_Norm)
          ##removed exp from columns
        } else {
          df <- df %>%
            select(ligand, key, construct, n, bottom, top, ec50, hill, ec50_sd, ec50_upper, ec50_lower, ec50_y_Norm)
        }
        
        
    
        #UPDATE SUBSET CONTROLS
        RV_df$ligands=unique(df$ligand) #Update Ligand List
        
        
          
        #SUBSET DF (REACTIVE)
        #Ligand
        if (!(input$ligand)=="") {
          sub_df <- df %>%
            filter(ligand==input$ligand)
          RV_df$ligkey = as.list(unique(sub_df$key)) #Update LigKey List
          # RV_df$constructs = as.list(unique(sub_df$construct)) #Update Construct List
        } else {
          sub_df <- df
          RV_df$ligkey=list()
          RV_df$constructs=list() #Update Construct List
        }
        
        
        #LigKey
        if (length(input$key>0)) {
          sub_df<-filter(sub_df, key %in% input$key)
          RV_df$constructs = as.list(unique(sub_df$construct)) #Update Construct List
        }
        
        #Construct
        if (length(input$construct>0)) {
          sub_df<-filter(sub_df, construct %in% input$construct)
          # RV_df$ligkey = as.list(unique(sub_df$key)) #Update LigKey List
        }
        
    
        
        #SUM_STAT SUBSET CONTROLS
        if (df_type=='auto') {
          #Plate_Type
          if (input$plate_type==2) {
            sub_df<-filter(sub_df, plate_type=='ensemble')
            RV_df$ligkey = as.list(unique(sub_df$key)) #Update LigKey List
            RV_df$constructs = as.list(unique(sub_df$construct)) #Update Construct List
          } else if (input$plate_type==3) {
            sub_df<-filter(sub_df, plate_type=='single')
            RV_df$ligkey = as.list(unique(sub_df$key)) #Update LigKey List
            RV_df$constructs = as.list(unique(sub_df$construct)) #Update Construct List
          }
        
          # Exp Range
          # RV_df$exp_range <- unique(sub_df$exp)
          # sub_df<-filter(sub_df, exp>=input$exp_range[1] & exp<=input$exp_range[2])
          
          
          #Success Rate Threshold
          sub_df <- filter(sub_df, success_rate>input$success_rate)
          
        }
        
        
        #PRODUCE FINAL PREVIEW DF
        RV_df$final_df <- sub_df
        

        
        
        
        
        ####DOSE RESPONSE CURVE
        if (nrow(sub_df)>1) {
        #Assign colors
        if (length(unique(df$construct)) <= 8) {
          construct.colors=brewer.pal(length(unique(df$construct)), "Dark2") #Dark2 works for 8 constructs in total, Set3 is one color larger
          names(construct.colors)<-levels(as.factor(unique(df$construct)))#Loads in constructs and assigns to color
        }
        
        #Create mini_DF to display population below graph
        pop_df<- as.data.frame(
          sub_df %>%
            select('Construct'=construct, 'n'=n, 'Key'=key, 'EC50'=ec50)
        )
        
        #Calculate Hill Equation (REACTIVE - Normalization & Subset)  
        hill_eq=data.frame()
        
        for(j in 1:nrow(sub_df)) {
          beta=as.numeric(sub_df[j,c('bottom', 'top', 'ec50', 'hill')])
    
          x_range=lseq(plot_bounds(sub_df)[1], plot_bounds(sub_df)[2], 500)
          
          if (input$normalize) {
            drug_curve=logistic_array(x = x_range, beta)/beta[2]*100
          } else {
            drug_curve=logistic_array(x = x_range, beta)
          }
          
          
          #Normalized to pcoda71
          # p71_top=eval(parse(text=i))[eval(parse(text=i))['construct']=='pcoda71', 'top']
          # drug_curve=logistic_array(x = x_range, beta)/p71_top
          
          construct=rep(sub_df[j,'construct'], 500)
          
          hill_eq=rbind(hill_eq, as.data.frame(cbind.data.frame(x_range, drug_curve, construct)))
          
          colnames(hill_eq)=c('x_range', 'drug_curve', 'construct')
        }
        
        #Plot Dose Response Curve  
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
          xlab(bquote(.(input$ligand) ~ "["*mu*M*"]")) +
          ylab('Normalized Peak Current \n (% of Imax)') +
          theme_grey(base_size = 22) +
          theme(legend.title = element_blank()) +
          ggtitle(paste(sub_df[j,2])) + #overwritten below
          ggtitle(label = NULL)
        
        # Subset rows that only have a standard dev greater than 1
        errorbar_df=sub_df[!is.na(sub_df['ec50_upper'] | sub_df['ec50_lower']),]
        
        drug_plot_stat= drug_plot_aes +
          geom_errorbarh(data=errorbar_df, mapping=aes(y=ec50_y_Norm, x=ec50, xmin=ec50_lower, xmax=ec50_upper, color=construct), height=7, size=1.5, alpha=0.5) +
          geom_point(data=sub_df, mapping=aes(group=construct, y=ec50_y_Norm, x=ec50, color=construct), size=5, shape=20, alpha=0.3, na.rm = TRUE)
        
        #Output DRC Plot
        if (input$geombars==TRUE) {
          RV_plot$drc=drug_plot_stat
        } else {
          RV_plot$drc=drug_plot_aes
        }
        
        #Output Mini_DF (n, Key. EC50 info)
        RV_plot$pop_df <- pop_df
        }
        
        
        
        
        ####Log10 FC Pipeline (REACTIVE)
        RV_FC$fc_constructs=sub_df$construct
        
        if (!(input$constructFC)=="") {
    
        #Select COnstruct for Normalization
        constructFC=input$constructFC
        constructFC_ec50=sub_df[sub_df$construct==constructFC, 'ec50']
        
        #Subset Dataframe
        FC_df <- sub_df %>% select(ligand, construct, n, 'log10FC_ec50'=ec50, 'log10FC_ec50sd'=ec50_sd)
    
        #log10 Transform and Normalize Dataframe
        FC_df$log10FC_ec50<-(log10(FC_df$log10FC_ec50)-log10(constructFC_ec50))/log10(constructFC_ec50) #Log10 transform and Normalization
        
        #Calculate Standard Deviation
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
        
        
        fc_gplot<-ggplot(data = FC_norm,
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
          ggtitle(paste0("Log10 Fold Change of EC50\n
                         (Normalized to ", constructFC, ")\n
                         [Ligand: ", input$ligand, "]"))+
          theme(plot.title = element_text(hjust = 0.5))
        
        #Update reactive value with FC Plot
        RV_FC$fc_plot=fc_gplot
      
    }
    
    }
    
  })

  #END OF PIPELINE
  
  
  
    
  #REACTIVE CONTROLS
  #Subset by Ligand
  observe({
    updateSelectInput(session, 'ligand', choices=RV_df$ligands, selected = "")
  })
  
  #Subset by Ligand Key
  observe({
    updateSelectizeInput(session, 'key', choices=RV_df$ligkey, selected = NULL)
  })
  
  #Subset by Construct
  observe({
    updateSelectizeInput(session, 'construct', choices=RV_df$constructs, selected = NULL,
                         options = list(
                           maxItems = 8, #Change this number if you want to adjust the total number of constructs you can subset by (will use default colorscheme if over 8)
                           placeholder = "Select or Type")
                         )
  })
  
  #subset by exp_date
  observe({
    updateDateRangeInput(session, 'exp_range', start = ymd(RV_df$exp_range[1]), end = ymd(tail(RV_df$exp_range, 1)))
  })

  
  # Select COnstruct to Normalize log10FC Plot to:
  observe({
    updateSelectInput(session, 'constructFC', choices=RV_FC$fc_constructs, selected = "")
  })

  
  #OUTPUT
  #Summary Statistics Dataframe
  output$final_df <- renderDT({
    datatable(RV_df$final_df, options = list(
      pageLength = 15,
      initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#2ad', 'color': '#fff'});",
        "}")
    ))
    
  })
  
  # Downloadable csv of selected dataset ----
  output$save_df <- downloadHandler(
    filename = function() {
      paste("gmean", input$ligand, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(RV_df$final_df, file, row.names = FALSE, col.names = T)
    }
  )
  
  #Dose Response Curve
  output$drc <- renderPlot({
    RV_plot$drc
  }, width = 750, height = 500)
  
  # Downloadable img of plot ----
  output$save_plot_png <- downloadHandler(
    filename = function() {
      paste0(input$ligand, "_drc.png")
    },
    content = function(file) {
      device <- function(..., width, height) grDevices::png(..., width = 10, height = 6.66, res = 300, units = "in")
      ggsave(file, plot=RV_plot$drc, device=device)
    }
  )
  

  #Mini DRC Dataframe (Gives Population and EC50 Values)
  output$pop_n <- renderTable({
    RV_plot$pop_df
  })
  
  
  #FC Plot
  output$fc_plot <- renderPlot({
    RV_FC$fc_plot
  }, width = 750, height = 500)
  
  # Downloadable img of plot ----
  output$save_FCplot <- downloadHandler(
    filename = function() {
      paste0(input$constructFC, "_normalized_log10FC.png")
    },
    content = function(file) {
      device <- function(..., width, height) grDevices::png(..., width = 10, height = 6.66, res = 300, units = "in")
      ggsave(file, plot=RV_plot$fc_plot, device=device)
    }
  )
  
  
})

#Written by Vasco Morais (March 2019) Cell: 415-845-2118