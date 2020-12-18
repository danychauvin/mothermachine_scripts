# 20201211, Dany Chauvin

# INSTALLING AND LOADING NECESSARY PACKAGES 

# Following packages are necessary to import data from deepMoma and analyze these data.

# WARNING: libstdc++.so.6 library, with proper GLIBCXX versions is necessary to use some packages imported here. 
# Be sure that GCCcore/8.3.0 is loaded/installed.
# To do so while using Rstudio on the scicore server "service06", do the following, add the following to your ~/.bashrc file:
# `if [[ "$HOSTNAME" = *service06* ]]; then
#     ml GCCcore/8.3.0
#  fi'

# Currently, packages are installed in my HOME folder, from centOS7 binaries at packagemanager.rstudio.com. I am not yet using renv, planned in the future.
# On Rstudio, set R version to 4.0.3.
# If this is the first time you use these codes, uncomment to following lines to install the corresponding packages.  

# Uncomment below when using for the first time.
options(repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/__linux__/centos7/latest"))
#install.packages("tidyverse")
#install.packages(c("here","cowplot"))
#install.packages("devtools")
#remotes::install_github(c('hadley/multidplyr'))
#install.packages("ccaPP")
#install.packages("~/R/vngMoM.tar", repos = NULL)
#remotes::install_github(c('julou/ggCustomTJ'))
#install.packages("reticulate")
#install.packages("ggcorrplot")
#install.packages("lemon")
#install.packages("LSD")
#install.packages("parallel")
#install.packages("ggpubr")
# Comment above if the packages are already installed.

# Necessary libraries
library(tidyverse)
library(RcppArmadillo)
library(tools)
library(here)
library(cowplot)
library(devtools)
library(multidplyr)
library(vngMoM)
library(ggCustomTJ)
library(renv)
library(svglite)
Sys.setenv(RETICULATE_PYTHON = '/scicore/home/nimwegen/rocasu25/.local/share/r-miniconda/envs/r-reticulate/bin/python')
library(reticulate)
library(ggcorrplot)
library(lemon)
library(parallel)
library(broom)
library(stats)
library(ggpubr)

#Set the conda environment
#Same as 'conda activate'
#use_condaenv('anaconda3') 

# Important functions that are used to import the data
# SET NECESSARY FUNCTIONS TO GENERATE PATHS TO PREPROCESSED DATA FILES ####
data2preproc_file <- function(.f)
  basename(.f) %>% sub("ExportedCellStats_", "", .) %>% 
  file_path_sans_ext %>% paste0("_frames.txt")
data2preproc <- function(.f)
  file.path(data2preproc_dir(.f), data2preproc_file(.f))


# FUNCTION TO SET OPTIMIZATION READY

set_optimization_ready <- function(data,cPath,fPath,fName,me,va,ste,num,N_cell){
  
  #Make the necessary folders
  dir.create(fPath,showWarnings=FALSE)
  dir.create(paste(fPath,"/","input",sep=""),showWarnings=FALSE)
  dir.create(paste(fPath,"/","output",sep=""),showWarnings=FALSE)
  #Copy paste the code: the code is always overwritten.
  system(paste("\\cp -r",cPath,fPath,sep=' '))

  # First, I will assume that in a given conditions, all growth lanes, all strains are alike.
  # Therefore, I will randomly sample 100 cell cycles.
  cell_sample_for_optimization <- data %>% 
    filter(medium==me) %>%
    group_by(cell) %>% 
    filter(n()>=num) %>% 
    ungroup() %>% 
    distinct(cell) %>% 
    sample_n(N_cell)
  
  sample_for_optimization <- data %>% 
    semi_join(cell_sample_for_optimization,by=c("cell"))
  
  #Adapting the format to the script
  sample_for_optimization <- sample_for_optimization %>%
    rename(cell_ori=cell) %>% 
    #parent_ID=parent_id) %>%
    mutate(pos=paste(pos,orientation,sep=".")) %>%  #NECESSARY CAUSE ATHOS CODE RECOMPUTES WRONGLY THE LANE_ID...
    mutate(pos=paste(pos,"_",sep="")) %>% #NECESSARY CAUSE ATHOS CODE RECOMPUTES WRONGLY THE LANE_ID...
    mutate(gl=paste(gl,"_",sep="")) %>% #NECESSARY CAUSE ATHOS CODE RECOMPUTES WRONGLY THE LANE_ID...
    mutate(cell=paste(lane_ID,id,sep="_"))
  #select(cell,date,time_sec,lane_ID,parent_ID,length_um)
  
  #Write sample
  readr::write_csv(sample_for_optimization,sprintf(paste(fPath,"/input/",fName,sep="")))
  
  #Append jobs to .input/commands.cmd
  lineToWrite <- sprintf("python %s/optimization_code/parameters_find.py %s/input/%s %s %s %s",fPath,fPath,fName,va,ste,num)
  system(sprintf('touch %s/input/commands.cmd',fPath))
  write(lineToWrite,file=sprintf('%s/input/commands.cmd',fPath),append=TRUE,N_cell)
  
  #Modify jobs.sh
  command_path <- sprintf("%s/input/commands.cmd",fPath) 
  results_path <- sprintf("%s/output/results_%%a",fPath)
  errors_path <- sprintf("%s/output/errors_%%a",fPath)
  jobs_path <- sprintf("%s/optimization_code/jobs.sh",fPath)
  system(sprintf("sed -i \"s+path_to_output+%s+g\" %s",results_path,jobs_path))
  system(sprintf("sed -i \"s+path_to_error+%s+g\" %s",errors_path,jobs_path))
  system(sprintf("sed -i \"s+path_to_commands+%s+g\" %s",command_path,jobs_path))
  
}

set_prediction_ready <- function(fPath,results_df,data){
  
  #Create necessary folders
  dir.create(paste(fPath,"/","prediction_input",sep=""),showWarnings=FALSE)
  dir.create(paste(fPath,"/","prediction_output",sep=""),showWarnings=FALSE)
  
  #Preparing the dataframe of conditions
  prediction_parameters <- results_df %>% 
    mutate(va="length_um") %>% 
    group_by(medium) %>% 
    arrange(log_likelihood) %>% 
    filter(row_number()==1) %>% #Here keeping only the parameters giving the best log likelihood.
    ungroup() %>% 
    arrange(medium) %>% 
    select(medium,va,ml,gamma,sl2,sm2,sd2) %>% 
    mutate(index=row_number()) %>% 
    mutate(input_path=sprintf("%s/prediction_input/%s_%s.csv",fPath,medium,as.character(index)),
           output_path=sprintf("%s/prediction_output/%s_%s.csv",fPath,medium,as.character(index)),
           lineToWrite=sprintf("python %s/optimization_code/path_predictions.py %s %s %s %s %s %s %s %s",fPath,input_path,output_path,va,ml,gamma,sl2,sm2,sd2))
  
  index_df <- prediction_parameters %>% 
    select(medium,index)
  
    #Copy paste data
    data %>%
      left_join(index_df,by=c("medium")) %>% 
      rename(cell_ori=cell) %>% 
      mutate(pos=paste(pos,orientation,sep=".")) %>%  #NECESSARY CAUSE ATHOS CODE RECOMPUTES WRONGLY THE LANE_ID...
      mutate(pos=paste(pos,"_",sep="")) %>% #NECESSARY CAUSE ATHOS CODE RECOMPUTES WRONGLY THE LANE_ID...
      mutate(gl=paste(gl,"_",sep="")) %>% #NECESSARY CAUSE ATHOS CODE RECOMPUTES WRONGLY THE LANE_ID...
      mutate(cell=paste(lane_ID,id,sep="_")) %>% 
      group_by(medium) %>% 
      do((function(.df){
        medium <- unique(.df$medium)
        readr::write_csv(.df,sprintf("%s/prediction_input/%s_%s.csv",fPath,medium,unique(.df$index)))
        return(data.frame())})(.))
  
  
  command_path <- sprintf("%s/prediction_input/commands.cmd",fPath)
  results_path <- sprintf("%s/prediction_output/results_%%a",fPath)
  errors_path <- sprintf("%s/prediction_output/errors_%%a",fPath)
  jobs_path <- sprintf("%s/optimization_code/prediction_jobs.sh",fPath)
  #system(sprintf('rm %s',command_path))
  system(sprintf('touch %s',command_path))
  
  #Adding jobs to command
  prediction_parameters %>% 
    arrange(index) %>% 
    group_by(index) %>% 
    do((function(.df){
      write(unique(.df$lineToWrite),file=sprintf('%s/prediction_input/commands.cmd',fPath),append=TRUE)
      return(data_frame())})(.)) %>% 
    ungroup()
  
  #Writing the prediction_jobs.sh
  command_path <- sprintf("%s/prediction_input/commands.cmd",fPath) 
  results_path <- sprintf("%s/prediction_output/results_%%a",fPath)
  errors_path <- sprintf("%s/prediction_output/errors_%%a",fPath)
  jobs_path <- sprintf("%s/optimization_code/prediction_jobs.sh",fPath)
  system(sprintf("sed -i \"s+path_to_output+%s+g\" %s",results_path,jobs_path))
  system(sprintf("sed -i \"s+path_to_error+%s+g\" %s",errors_path,jobs_path))
  system(sprintf("sed -i \"s+path_to_commands+%s+g\" %s",command_path,jobs_path))
  
  runAllJobsLine <- sprintf("sbatch --array=1-$(cat %s|wc -l):1 prediction_jobs.sh",command_path)
  print(runAllJobsLine)
  return(prediction_parameters)
}


concatenate_traces <- function(df){
  # Input: should have a growth_rate field, depending on the one that you want to consider for the computation.
  # Output: concatenated cell cycles
  cell_and_parent_1 <- df %>% 
    distinct(cell,.keep_all=TRUE) %>% 
    select(cell,parent,medium) %>% 
    rename(degree_1=parent)
  
  cell_and_parent_2 <- df %>%
    distinct(cell,.keep_all = TRUE) %>% 
    select(cell,parent,medium) %>% 
    rename(degree_1=cell) %>% 
    rename(degree_2=parent)
  
  cell_and_parent_3 <- df %>%
    distinct(cell,.keep_all = TRUE) %>% 
    select(cell,parent,medium) %>% 
    rename(degree_2=cell) %>% 
    rename(degree_3=parent)
  
  cell_and_parent_4 <- df %>%
    distinct(cell,.keep_all = TRUE) %>% 
    select(cell,parent,medium) %>% 
    rename(degree_3=cell) %>% 
    rename(degree_4=parent)
  
  phylogenies <- left_join(cell_and_parent_1,cell_and_parent_2,by=c("degree_1","medium")) %>% 
    left_join(cell_and_parent_3,by=c("degree_2","medium")) %>% 
    left_join(cell_and_parent_4,by=c("degree_3","medium")) %>% 
    mutate(degree_0=cell) %>%
    rename(trace_ref=cell) %>% 
    select(trace_ref,degree_0,degree_1,degree_2,degree_3,degree_4,medium) %>% #here limited to 3 consecutive growth
    gather(degree,cell_id,c(degree_0,degree_1,degree_2,degree_3,degree_4)) %>%  
    mutate(degree=as.double(substr(degree,nchar(degree),nchar(degree)))) %>% 
    arrange(trace_ref,desc(degree)) %>% 
    rename(cell=cell_id)
  
  # Now join phylogenies to another df
  concatenatedTraces <- phylogenies %>% 
    left_join(df,by=c("cell","medium")) %>% 
    drop_na() %>% #Dropping cells that do not exist in the dataset
    select(trace_ref,degree,cell,time_sec,medium,growth_rate) %>% 
    arrange(trace_ref,time_sec)
  
  return(concatenatedTraces)
}

autocorrelationFunction <- function(full_df,mean_div_time,dt,nCellCycle,cond){
  results <- c()
  M <- round(mean_div_time/dt*nCellCycle,0) #Autocorrelation is computed over nCellCycle * mean number of data point per cell cycle
  for (i in 1:M){
    results[i] <- autocorrelationFunction_lag(full_df,i)}
  results_df <- tibble(condition=cond,lag_step=c(0:(M-1)),lag_min=c(0:(M-1))*dt,autocorrelation=results)
  return(results_df)
}
  
autocorrelationFunction_lag <- function(df,lag){
  #input: df with trace_ref,cell,time_sec and growth_rate, lag as an integer
  #output: autocorrelation value for a given lag.
  
  new_df <- df %>% 
    mutate(uid=paste(cell,time_sec,sep=".")) %>% 
    group_by(trace_ref) %>% 
    arrange(time_sec) %>%
    select(trace_ref,uid,growth_rate) %>% 
    do((function(.df){
      N <- nrow(.df)
      if (lag>N){
      return(tibble())
        }else{
      concat <- cbind(
        .df[c(lag:N),],
        .df[c(1:(N-lag+1)),])
    colnames(concat) <- list('trace_ref','uid_1','growth_rate_1','trace_ref_2','uid_2','growth_rate_2')
    return(concat)}
    })(.)) %>% 
    ungroup() %>%
    select(-c(trace_ref,trace_ref_2)) %>% 
    unique()
  corvalue <- (mean(new_df$growth_rate_1*new_df$growth_rate_2)-mean(new_df$growth_rate_1)*mean(new_df$growth_rate_2))/(sd(new_df$growth_rate_1)*sd(new_df$growth_rate_2))
  return(corvalue)}

myframes_mother_daughters_compute <- function(df){
  #input:
  #cell: with the format paste(date,pos,orientation,gl,sep=".")
  #time_sec: time in sec since the beginning of the experiment
  #condition: before the switch it should be "glucose_pre_switch"
  #description: should match the medium after the switch
  #parent_id: should be a integer
  #id: shoud be an integer
  #length_um
  #gfp_nb
  
  #output:
  #description: condition experienced by the cells after the switch
  #condition: condition experienced by the cells at time_sec
  #trace_id: 'cell' of the mother cell that experience the switch
  #type: mother or daughter
  #time_sec: idem above
  #length_um: if mother, simply the gfp_nb of the mother, if daughters: sum of length_um
  #gfp_nb: if mother, simply the gfp_nb of the mother, if daughters: sum of gfp_nb
  
  cell_experiencing_switch <- df %>% 
  group_by(cell) %>% 
  arrange(time_sec) %>% 
  filter(row_number() %in% c(1,n())) %>% #Cells that have frames in "glucose pre-shift and condition==description are kept"
  do((function(.df){
    test <- (unique(.df[1,]$condition)=="glucose_pre_switch")*(unique(.df[2,]$condition) %in% c(unique(.df[2,]$description),"glucose_post_switch"))
    if(test){
      return(.df[1,])}
    else{
      return(tibble())}})(.)) %>% 
  ungroup() %>% 
  mutate(parent=cell) %>% 
  select(cell,parent)

myframes_switch <- df %>%  #cells that experience the switch
  semi_join(cell_experiencing_switch,by=c("cell")) %>% 
  mutate(trace_id=cell) %>% 
  mutate(type="mother") %>% 
  ungroup() %>% 
  select(description,condition,trace_id,type,time_sec,length_um,gfp_nb)
  
daughters <- df %>% 
  mutate(parent=paste(date,pos,orientation,gl,parent_id,sep=".")) %>% 
  semi_join(cell_experiencing_switch,by=c("parent")) %>% 
  group_by(parent) %>%
  arrange(cell) %>% 
  mutate(daughter_1=ifelse(row_number()==1,cell,NA)) %>% 
  fill(daughter_1,.direction=c("down")) %>% 
  mutate(daughter_num=ifelse(cell==daughter_1,1,2)) %>% 
  select(-daughter_1) %>% 
  ungroup() %>% 
  select(description,condition,parent,daughter_num,time_sec,length_um,gfp_nb) %>% 
  group_by(parent) %>% 
  do((function(.df){
    mat1 <- as.matrix(.df[c("length_um","gfp_nb","daughter_num")] %>% filter(daughter_num==1) %>% select(-daughter_num))
    mat2 <- as.matrix(.df[c("length_um","gfp_nb","daughter_num")] %>% filter(daughter_num==2) %>% select(-daughter_num))
    n <- min(nrow(mat1),nrow(mat2))
    if(n==0){return(tibble())}
    else{
      #print(n)
      mat <- mat1[c(1:n),]+mat2[c(1:n),]
      if(n==1){mat <- as_tibble_row(mat)} #bad as_tibble behaviour...
      else{mat <-as_tibble(mat)}
      #print(mat)
      infos <- .df %>% 
        filter(daughter_num==1) %>% 
        select(description,condition,parent,time_sec)
      #print(colnames(infos))
      new_df <- cbind(infos[c(1:n),],mat) #if nrow(mat)<2, as_tibble messes up the format: columns no longer exists, only value stays...
      #if(length(colnames(new_df))!=6){print(unique(new_df$parent))}
      colnames(new_df) <- list("description","condition","parent","time_sec","length_um","gfp_nb")
      return(new_df)}})(.)) %>% 
  ungroup() %>% 
  rename(trace_id=parent) %>% 
  mutate(type="daughters")
  
new_df <- rbind(
  myframes_switch,daughters) %>% 
  arrange(trace_id,time_sec) %>% 
  select(description,condition,trace_id,type,time_sec,length_um,gfp_nb)

return(new_df)

}
