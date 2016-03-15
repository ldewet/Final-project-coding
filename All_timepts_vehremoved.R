together_timepoints <- function(datafile, templatefile, plottitle, graphname){
  
  #import dplyr
  library(dplyr)
  
  #Set working directory to folder with files
  #setwd("/Users/Lari/Documents/UChicago/Winter 2016/Computing")
  
  #read data
  PCRdata <- read.table("Data.csv", sep = ",", na.strings = "N/A", stringsAsFactors = FALSE)
  template <- read.table("Template_1V4.csv", stringsAsFactors = FALSE, header = TRUE, sep = ",", 
                         row.names = 1)
  
  #delete the headers
  PCRdata <- PCRdata[-1,]
  #change column names to numbers
  colnames(template) <- c(1,2,3,4,5,6,7,8,9,10,11,12)
  
  #now let's split the wells into number and letter parts
  lett_vec <- vector(mode = 'character', length = nrow(PCRdata))
  num_vec <- vector(mode = 'numeric', length = nrow(PCRdata))
  
  #Now assign each of the letters from our data into the vector
  for(i in 1:nrow(PCRdata)){
    lett_vec[i] <- strsplit(PCRdata$V1[i], split = '[0-9]{1,2}')
  }
  lett_vec <- unlist(lett_vec)
  
  
  #And assign the numbers from our data into the number vector
  for(i in 1:nrow(PCRdata)){
    #num_vec[i] <- strsplit(PCRdata$V1[i], split = '[A-Z]')
    num_vec[i] <- as.numeric(gsub("[^\\d]+", "", PCRdata$V1[i], perl=TRUE))
  }
  
  #add these vectors to the dataframe
  PCRdata$letter <- lett_vec
  PCRdata$number <- num_vec
  
  #rename the headers
  names(PCRdata) <- c('Well', 'Threshold_cycle', 'letter', 'number')
  
  #Add the gene name, treatment, timepoint and whether it's a housekeeping gene
  genes <- vector(mode = 'character', length = nrow(PCRdata))
  treatment <- vector(mode = 'character', length = nrow(PCRdata))
  timept <- vector(mode = 'character', length = nrow(PCRdata))
  housekeeping <- vector(mode = 'character', length = nrow(PCRdata))
  PCRdata$gene <- genes
  PCRdata$treatment <- treatment
  PCRdata$timept <- timept
  PCRdata$HK <- housekeeping
  
  #now add the strings in. If housekeeping gene, value = 1
  for(i in 1:nrow(PCRdata)){
    rowTF <- row.names(template) == PCRdata$letter[i]
    fromtemplate <- dplyr::select(template, PCRdata$number[i])  %>% dplyr::filter(rowTF)
    splitstrings <- strsplit(as.character(fromtemplate), "-")[[1]]
    PCRdata$gene[i] <- splitstrings[1]
    PCRdata$treatment[i] <- splitstrings[2]
    PCRdata$timept[i] <- splitstrings[3]
    PCRdata$HK[i] <- splitstrings[4]
  }
  
  #convert threshold cycle values to numeric (currently is character- dplyr::summarise can't use this)
  PCRdata$Threshold_cycle <- as.numeric(PCRdata$Threshold_cycle)
  
  #turn into character vectors
  PCRdata$treatment <- as.character(PCRdata$treatment)
  #turn the treatments into factors so they stay in order in the legend
  PCRdata$treatment <- factor(PCRdata$treatment, levels = PCRdata$treatment)
  
  #group by gene, treatment and timepoint for all non-housekeeping genes
  summary1 <- dplyr::filter(PCRdata, HK == 0) %>% dplyr::group_by(gene, treatment, timept) %>% na.omit %>% dplyr::summarise(avg = mean(Threshold_cycle))
  #group by gene, treatment and timepoint for all housekeeping genes
  summary2 <- dplyr::filter(PCRdata, HK == 1) %>% dplyr::group_by(gene, treatment, timept) %>% na.omit %>% dplyr::summarise(HK_avg = mean(Threshold_cycle))
  #now let's merge these together
  summary_join <- dplyr::full_join(summary1, summary2)
  
  #find the unique timepts
  uniquetimepts <- vector(mode = "character")
  uniquetimepts <- unique(PCRdata$timept)
  uniquetimepts <- uniquetimepts[!is.na(uniquetimepts)]
  
  #let's split by timepoint
  joinedlist <- list()
  for(i in uniquetimepts){
    #assign each timept to a separate dataframe
    assign(i, data.frame()) 
    #filter each timept and the values to the datafram
    df_tmpt <- dplyr::filter(summary_join, timept == i) 
    #normalise the gene to the housekeeping gene
    sub_HK <- dplyr::mutate(df_tmpt, normalisedgene = avg - HK_avg)
    #subract the vehicle control values from each of the treatments
    sub_veh <- dplyr::mutate(sub_HK, subtractveh = normalisedgene - sub_HK$normalisedgene[1])
    #calculate the fold change
    fchange <- dplyr::mutate(sub_veh, foldchange = (2 ^ -(subtractveh)))
    #put all the dataframes in a list so that we can merge them
    joinedlist[[i]] <- fchange
  }
  
  #let's merge these dataframes back together so we can graph with all timepts
  joined_df <- do.call("rbind", joinedlist)
  
  veh_removed <- dplyr::filter(joined_df, treatment != joined_df$treatment[1]) 
  
  #make a pretty graph
  library(ggplot2)
  
  #turn into character vectors
  veh_removed$treatment <- as.character(veh_removed$treatment)
  #turn the treatments/timepoints into factors so they stay in order in the legend
  veh_removed$treatment <- factor(veh_removed$treatment, levels = veh_removed$treatment)
  veh_removed$timept <- factor(veh_removed$timept, levels = veh_removed$timept)
  
  
  #finally, let's plot it!!
  graph <- ggplot(data=veh_removed, aes(x=interaction(timept), y=foldchange, fill = treatment)) + 
    geom_bar(colour="black", aes(fill = treatment), stat="identity", position = position_dodge()) + 
    ggtitle(plottitle) + 
    theme(plot.title = element_text(face="bold", size=20), axis.title.x = element_text(size = 16, face = "bold"), 
          axis.text.x = element_text(size = 15), axis.title.y = element_text(size = 14, angle=90, face = "bold"), 
          axis.text.y = element_text(size = 13), panel.background = element_rect(fill = "white",
                                                                                 colour = "white", size = 0.5, linetype = "solid"), panel.border = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(size = 0.5, linetype = "solid", colour = "black")) + 
    labs(y= "Fold Change over Vehicle", x = "Treatment") +
    scale_fill_manual(values = c("#66CC00", "#000066", "#66CCCC", "#990099", "black"), name="Treatment") 
  
  #save the graph to a jpg file in the working directory
  ggsave(paste(graphname, "jpg", sep = "."))
}

