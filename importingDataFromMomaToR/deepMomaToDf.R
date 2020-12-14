# IMPORT DEEP MOMA DATA INTO R DATAFRAME
# Author: Dany Chauvin
# Date: 20201125

# Input
# Path to the dataframe in the .csv file which recapitulates experimental details.
# Ouput
# mycells and myframes_to_mycells dataframes

# SET DEFAULT VARIABLES
dl <- 0.065 # µm/pixel
vertical_cutoff <- 150 #Necessary to get rid of cells that are touching the top of the growth lane where segmentation is not that good anymore.

# Parallel environment for multidplyr
mycluster <- min(30, parallel::detectCores()-1) %>%  # do not use more than 30 cores
  default_cluster() %>% cluster_library( # load currently loaded packages on each core
    names(sessionInfo()$otherPkgs))

myconditions <- readr::read_csv(path_to_data_summary,
                               col_types = cols(
                                 date=col_character(),
                                 description=col_character(),
                                 f_start=col_double(),
                                 f_end=col_double(),
                                 condition=col_character(),
                                 t_interval=col_double(),
                                 orientation=col_character(),
                                 data_path=col_character()))
                                                

# Computing t_start and t_end in seconds
myconditions <- myconditions %>% 
  group_by(date,description,data_path) %>% 
  arrange(f_start) %>%
  mutate(t_end=ifelse(row_number()==1,(f_end)*t_interval*60,NA)) %>% 
  mutate(t_start=ifelse(row_number()==1,(f_start)*t_interval*60,NA)) %>% 
  do((function(.df){
    new_df <- .df
    if(nrow(new_df)>=2){
      for(i in c(2:nrow(new_df))) 
      {
        new_df$t_start[i]=new_df$t_end[i-1]
        new_df$t_end[i]=new_df$t_start[i]+(new_df$f_end[i]-new_df$f_start[i])*new_df$t_interval[i]*60}
    }
    return(new_df)
  })(.)) %>%
  ungroup()

myconditions <- myconditions %>% 
  group_by(date,description,data_path) %>% 
  do((function(.df){
    new_df <- find.files(unique(.df$data_path), "ExportedCellStats_*.csv") %>% 
      data.frame(file_path=., stringsAsFactors=FALSE)
    new_df <- crossing(.df,new_df)
    return(new_df)
  } )(.)) %>% 
  ungroup()

nFiles <- myconditions %>% distinct(file_path) %>% nrow()

#nc <- min(nFiles, length(mycluster)) # this is a dirty hack because multidplyr crashes with less shards than cores
#Proper column names
properColNames <- c("file_path","lane_ID","cell_ID","frame","cell_rank","genealogy","type_of_end","parent_ID","cells_in_lane","bbox_top px","bbox_bottom px","center_x px","center_y px","width px","length px","tilt rad","area px^2","bgmask_area px^2","fluo_cellmask_1","fluo_bgmask_ch_1","fluo_ampl_ch_1","fluo_bg_ch_1")

myframes <- myconditions %>% 
  ungroup() %>% 
  distinct(file_path) %>% 
  .$file_path %>% 
  lapply(function(.l) readr::read_delim(.l,skip=2,delim=";") %>% 
           mutate(file_path=.l)) %>% 
  lapply(function(.l) .l %>% select(all_of(properColNames))) %>% 
           do.call(rbind, .)
         
myframes <- myframes %>%
  left_join(myconditions,by=c('file_path')) %>% 
  group_by(file_path,condition) %>% 
  filter(between(frame,f_start,f_end-1)) %>% 
  extract(file_path, c("date","pos", "gl"), ".*(\\d{8})_.*[Pp]os(\\d+).*_GL(\\d+).*", remove=FALSE, convert=TRUE) %>%
  ungroup()

myframes <- myframes %>% 
  rename(id=cell_ID,
         end_type=type_of_end,
         parent_id=parent_ID,
         vertical_bottom='bbox_bottom px',
         vertical_top='bbox_top px',
         center_x_pixel='center_x px',
         center_y_pixel='center_y px',
         width_pixel='width px',
         length_pixel='length px',
         tilt_radian='tilt rad',
         area_pixel2='area px^2',
         background_mask_aread_pixel2='bgmask_area px^2',
         fluo_cell_mask_ch1='fluo_cellmask_1',
         fluo_background_mask_ch1='fluo_bgmask_ch_1',
         fluo_amplitude='fluo_ampl_ch_1',
         fluo_background_ch1='fluo_bg_ch_1') %>%
  #Convert date to factors
  mutate_at(vars(date), factor) %>% 
  #Set cell ref
  mutate(cell=paste(date,pos,orientation,gl,id,sep='.'),
         parent=paste(date,pos,orientation,gl,parent_id,sep='.'),
         lane_ID=paste(date,pos,orientation,gl,sep='.')) %>% 
  #Compute length in um
  mutate(length_um=length_pixel*dl,
         width_um=width_pixel*dl) %>%
  #Compute vertical center
  mutate(vertical_center=(vertical_bottom + vertical_top)/2) %>% 
  # Propagating end_type information
  mutate(gl_id=paste(date, pos, orientation , gl, sep='.')) %>% 
  group_by(cell) %>% 
  fill(end_type,.direction="up") %>%
  mutate(discard_top=(vertical_top < vertical_cutoff)) %>%
  mutate(end_type=ifelse(any(discard_top), 'touchtop', end_type)) %>%
  ungroup()

myframes <- myframes %>% 
  mutate(time_sec=t_start+(frame-f_start)*t_interval*60) %>% 
  group_by(cell) %>% 
  partition(cluster=mycluster) %>%
  mutate(start_time=first(time_sec), end_time=last(time_sec),
         #cell_cycle=ifelse(end_type=='div', ((1:length(frame))-.5) / length(frame), NA),
         length_predict=fit_exp_elongation(time_sec, length_um)) %>% #using Michael ellipse fitting
  collect() %>% 
  arrange(cell, frame) %>% # sort data after `partition()`
  ungroup()

# Add the mean_position of the cell between two divisions in the growth lane.
myframes <- myframes %>% 
  group_by(cell) %>%
  mutate(mean_position_um=mean(vertical_center)*dl) %>% 
  ungroup()

# COMPUTE INSTANTANEOUS RATES ####
# discrete fluo derivative
#myframes <- myframes %>% 
#  group_by(date, pos, orientation, gl, id) %>% 
#  mutate(gfp_deriv=(gfp_nb-lag(gfp_nb))/(time_sec-lag(time_sec)),
#         l_deriv=(length_um-lag(length_um))/(time_sec-lag(time_sec))) %>% 
#  ungroup

# Add column switch to cells that are experiencing a switch
## LARGE NUMBER OF CELLS LOST HERE
myframes <- myframes %>%
  group_by(cell) %>% 
  arrange(frame) %>% 
  mutate(switch=ifelse(any((time_sec-t_start)<2*3600),TRUE,FALSE)) %>% 
  ungroup()

# Compute mycells only for...
# Cell whose parents are known, that divides, that are characterized by more than 4 data points
# Possible to discard cells that also experienced a switch at this point

myframes_to_mycells <- myframes %>% 
  group_by(cell) %>% 
  filter(parent_id!=-1) %>% 
  filter(end_type=='div') %>% 
  #filter(switch==FALSE) %>% 
  filter(n()>4) %>% 
  ungroup()

mycells <- myframes_to_mycells %>% 
  group_by(date,pos,orientation,gl,cell,description,condition) %>% 
  partition(cluster=mycluster) %>%
  do((function(.df) {
    .mod_ll_t <- fastLmPure( cbind(1, .df$time_sec), log(.df$length_um) )
    .mod_l_t <- fastLmPure( cbind(1, .df$time_sec), .df$length_um )
    .time_birth <- first(.df$time_sec)
    .time_div <- last(.df$time_sec)
    data.frame(npoints=.mod_ll_t$df.residual+1,
               time_birth=.time_birth, time_div=.time_div, div_time=.time_div-.time_birth,
               l_birth=first(.df$length_predict), l_div=last(.df$length_predict),
               logl_time_slope=.mod_ll_t$coefficients[2], logl_time_slopesd=.mod_ll_t$stderr[2],
               logl_time_r2=cor(.df$time_sec, log(.df$length_um))^2,
               l_time_slope=.mod_l_t$coefficients[2], l_time_slopesd=.mod_l_t$stderr[2], 
               l_time_r2=cor(.df$time_sec, .df$length_um)^2)
  })(.) ) %>% 
  collect() %>% 
  #arrange(condition, date, pos, gl, id) %>% 
  ungroup()

mycells <- mycells %>% 
  # filter by r2 of exponential fit
  filter(logl_time_r2>0.95) %>% 
  # create new variables of interest
  mutate(dl=l_div - l_birth,
         alpha=log(l_div/l_birth) / div_time) # ok since l_xxx are fitted
         #c_birth=g_birth/l_birth,
         #c_div=g_div/l_div,
         #dg=g_div - g_birth,
         #dcdt=(g_div/l_div - g_birth/l_birth) / div_time,
         #g=log(g_div/g_birth) / div_time, # would be cleaner from a fit but difficult to compare
         #gamma=dg/dl,
         #q=dg/dl*alpha)
  
myframes %>% distinct(cell) %>% nrow()
myframes_to_mycells %>% distinct(cell) %>% nrow()
mycells %>%  distinct(cell) %>% nrow()
myframes_to_mycells <- myframes_to_mycells %>% 
  semi_join(mycells,by=c("cell"))



           