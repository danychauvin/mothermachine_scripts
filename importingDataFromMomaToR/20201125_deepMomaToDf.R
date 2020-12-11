# DOCUMENTATION
# Author: Dany Chauvin
# Date: 20201125

# TODO/Improvment
# Set the vertical_cutoff right for deepMoma #vertical_cutoff <- 130#4 / dl # pixels
# In the future, the orientation should be added to the .csv file name during export by deep moma.
# duration is misleading: this is in fact switch_frame to change everywhere once finished
# medium: should be condition
# condition: change to "experiment_type"

# SET DEFAULT VARIABLES
dl <- 0.065 # µm/pixel
vertical_cutoff <- 150 #All cells characterized by shrinking end of cell cycle are there!
dir.create(here("slogs"), showWarnings=FALSE) # create a directory to store logs from the queue

# set a parallel environment to run multidplyr
mycluster <- min(30, parallel::detectCores()-1) %>%  # do not use more than 30 cores
  default_cluster() %>% cluster_library( # load currently loaded packages on each core
    names(sessionInfo()$otherPkgs))

# myconditions: dataframe of experimental conditions ####
# Each experiment is defined by, path, date,condition,duration,medium,t_interval, orientation in the dataframe: my conditions.
# path: should be an absolute path that leads to the deepMoma curated .csv files
# e.g: "/scicore/home/nimwegen/GROUP/MM_Data/Dany/20190730/202008_deep/20190730_chr_rpsB_curated_to_R"
# The folder name includes the strain (here promoter reporter strain "rpsB"), the vector ("chr"). Which are extracted by the script and added to the final dataframe.
# All the .csv files in that folder should have the same date, condition, duration, medium, t_interval and orientation variables.
# The .csv files should have a name such as: "ExportedCellStats_20190730_rpsB_rrnB_plac_6300_glu_gly_4_MMStack_Pos27_preproc_GL03_curated.csv"
# The name of the .csv file includes the position ("Pos27") and the growth lane ("GL03"). These are extracted by the script and added to the final dataframe.
# date: character format such as: "20201125"
# condition: character description of the type of experiment, e.g: "glu_lac_switch"
# duration: vector containing frames at which conditions change in the experiment (medium, t_interval, etc.). The last condition "lasts forever" and is therefore Inf. The length of this vector should be the same as medium and t_interval.
# medium: vector containing successive media encountered by the cells in the experiment.
# t_interval: vector that contains successive acquisition time periods during the experiment.

myconditions <- rbind(
  
  tibble(date='20200728',
         condition='glu_aa',
         f_start=0,
         f_end=Inf,
         medium='glu_aa',
         t_interval=1.5,
         orientation='b') %>% 
    crossing(tibble(data_path=c("./data/Dany/20200728/20200728_bottom_chr_hi1_curated",
                                "./data/Dany/20200728/20200728_bottom_chr_hi3_curated",
                                "./data/Dany/20200728/20200728_bottom_chr_med3_curated",
                                "./data/Dany/20200728/20200728_bottom_chr_rpsB_curated"))),
  
  tibble(date='20200728',
         condition='glu_aa',
         f_start=0,
         f_end=Inf,
         medium='glu_aa',
         t_interval=1.5,
         orientation='t') %>% 
    crossing(tibble(data_path=c("./data/Dany/20200728/20200728_top_chr_hi1_curated",
                                "./data/Dany/20200728/20200728_top_chr_rpsB_curated"))),
  
  tibble(date='20200730',
         condition='glu_aa',
         f_start=0,
         f_end=Inf,
         medium='glu_aa',
         t_interval=1.5,
         orientation='b') %>% 
    crossing(tibble(data_path=c("./data/Dany/20200730/20200730_bottom_chr_rplN_curated",
                                "./data/Dany/20200730/20200730_bottom_chr_rpmB_curated",
                                "./data/Dany/20200730/20200730_bottom_chr_rrnB_curated",
                                "./data/Dany/20200730/20200730_bottom_chr_med2_curated"))),
  
  tibble(date='20200730',
         condition='glu_aa',
         f_start=0,
         f_end=Inf,
         medium='glu_aa',
         t_interval=1.5,
         orientation='t') %>% 
    crossing(tibble(data_path=c("./data/Dany/20200730/20200730_top_chr_rplN_curated",
                                "./data/Dany/20200730/20200730_top_chr_rpmB_curated",
                                "./data/Dany/20200730/20200730_top_chr_rrnB_curated",
                                "./data/Dany/20200730/20200730_top_chr_med2_curated"))),
  
  tibble(date='20190515',
         condition='glu_gly_switch',
         f_start=c(0,440),
         f_end=c(440,Inf),
         medium=c('glucose','glycerol'),
         t_interval=c(3,6),
         orientation='b') %>% 
    crossing(tibble(data_path=c("./data/Dany/20190515/202008_deep/20190515_chr_hi1_curated_to_R",
                                "./data/Dany/20190515/202008_deep/20190515_chr_med2_curated_to_R",
                                "./data/Dany/20190515/202008_deep/20190515_chr_rpmB_curated_to_R"))),
  
  tibble(date='20190730',
         condition='glu_gly_switch',
         f_start=c(0,240),
         f_end=c(240,Inf),
         medium=c('glucose','glycerol'),
         t_interval=c(3,6),
         orientation='b') %>% 
    crossing(tibble(data_path=c("./data/Dany/20190730/202008_deep/20190730_chr_rpsB_curated_to_R",
                                "./data/Dany/20190730/202008_deep/20190730_chr_rrnB_curated_to_R"))),
  
  tibble(date='20190529',
         condition='glu_gly_switch',
         f_start=c(0,220),
         f_end=c(220,Inf),
         medium=c('glucose','glycerol'),
         t_interval=c(3,6),orientation='b') %>% 
    crossing(tibble(data_path=c("./data/Dany/20190529/202008_deep/20190529_chr_rplN_curated_to_R",
                                "./data/Dany/20190529/202008_deep/20190529_chr_med3_curated_to_R",
                                "./data/Dany/20190529/202008_deep/20190529_chr_hi3_curated_to_R"))),
  
  tibble(date='20200812',
         condition='ace',
         f_start=c(0),
         f_end=c(116),
         medium=c('ace'),
         t_interval=c(18.75),
         orientation='b') %>% 
    crossing(tibble(data_path=c("./data/Dany/20200812/20200812_bottom_chr_rpsB_curated",
                                "./data/Dany/20200812/20200812_bottom_chr_rplN_curated",
                                "./data/Dany/20200812/20200812_bottom_chr_rpmB_curated",
                                "./data/Dany/20200812/20200812_bottom_chr_rrnB_curated",
                                "./data/Dany/20200812/20200812_bottom_chr_med3_curated",
                                "./data/Dany/20200812/20200812_bottom_chr_hi1_curated",
                                "./data/Dany/20200812/20200812_bottom_chr_hi3_curated"))),
  
  tibble(date='20200812',
         condition='ace',
         f_start=c(0),
         f_end=c(116),
         medium=c('ace'),
         t_interval=c(18.75),
         orientation='t') %>% 
    crossing(tibble(data_path=c("./data/Dany/20200812/20200812_top_chr_rpsB_curated",
                                "./data/Dany/20200812/20200812_top_chr_rrnB_curated",
                                "./data/Dany/20200812/20200812_top_chr_med3_curated",
                                "./data/Dany/20200812/20200812_top_chr_hi1_curated",
                                "./data/Dany/20200812/20200812_top_chr_hi3_curated"))),
  
  tibble(date='20200911',
         condition='ace',
         f_start=c(0),
         f_end=c(Inf),
         medium=c('ace'),
         t_interval=c(18.75),
         orientation='b') %>% 
    crossing(tibble(data_path=c("./data/Dany/20200911/20200911_bottom_chr_med2_curated"))),
  
  tibble(date='20200911',
         condition='ace',
         f_start=c(0),
         f_end=c(Inf),
         medium=c('ace'),
         t_interval=c(18.75),
         orientation='t') %>% 
    crossing(tibble(data_path=c("./data/Dany/20200911/20200911_top_chr_med2_curated")))
  
) #End of rbind

# Computing t_start and t_end in seconds
myconditions <- myconditions %>% 
  group_by(date,condition,data_path) %>% 
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
  group_by(date,condition,data_path) %>% 
  do((function(.df){
    new_df <- find.files(unique(.df$data_path), "ExportedCellStats_*.csv") %>% 
      data.frame(file_path=., stringsAsFactors=FALSE)
    new_df <- crossing(.df,new_df)
    return(new_df)
  } )(.)) %>% 
  ungroup()

nFiles <- myconditions %>% distinct(file_path) %>% nrow()

nc <- min(nFiles, length(mycluster)) # this is a dirty hack because multidplyr crashes with less shards than cores
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
  group_by(file_path,medium) %>% 
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
  partition(cluster=mycluster[1:nc]) %>%
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

# Sanity check
#myframes %>% filter(date=="20190730") %>% filter(frame==240) %>% .$cell
#myframes %>% filter(cell=='20190730.28.b.9.91') %>% select(cell,medium,frame,time_sec,end_type) %>% View() 

myframes <- myframes %>% 
  ungroup() %>% 
  left_join(bind_rows(
   
    myframes %>% ungroup %>% filter(date=="20200728",orientation=="b") %>% 
      group_by(gl_id,condition) %>% summarise(file_path=unique(file_path)) %>% 
      extract(file_path, "strain", "./data/Dany/20200728/20200728_bottom_(.*)_curated/20200728*"),
    
    myframes %>% ungroup %>% filter(date=="20200728",orientation=="t") %>% 
      group_by(gl_id,condition) %>% summarise(file_path=unique(file_path)) %>% 
      extract(file_path, "strain", "./data/Dany/20200728/20200728_top_(.*)_curated/20200728*"),
    
    myframes %>% ungroup %>% filter(date=="20200730",orientation=="b") %>% 
      group_by(gl_id,condition) %>% summarise(file_path=unique(file_path)) %>% 
      extract(file_path, "strain", "./data/Dany/20200730/20200730_bottom_(.*)_curated/20200730*"),
    
    myframes %>% ungroup %>% filter(date=="20200730",orientation=="t") %>% 
      group_by(gl_id,condition) %>% summarise(file_path=unique(file_path)) %>% 
      extract(file_path, "strain", "./data/Dany/20200730/20200730_top_(.*)_curated/20200730*"),
    
    # 20190515: mxMoM (separated folders per strain)
    myframes %>% ungroup %>% filter(date=="20190515") %>% 
      group_by(gl_id,condition) %>% summarise(file_path=unique(file_path)) %>% 
      extract(file_path, "strain", "./data/Dany/20190515/202008_deep/20190515_(.*)_curated_to_R*"),
    
    # 20190529: mxMoM (separated folders per strain)
    myframes %>% ungroup %>% filter(date=="20190529") %>% 
      group_by(gl_id,condition) %>% summarise(file_path=unique(file_path)) %>% 
      extract(file_path, "strain", "./data/Dany/20190529/202008_deep/20190529_(.*)_curated_to_R*"),
    
    # 20190730: mxMoM (separated folders per strain)
    myframes %>% ungroup %>% filter(date=="20190730") %>% 
      group_by(gl_id,condition) %>% summarise(file_path=unique(file_path)) %>% 
      extract(file_path, "strain", "./data/Dany/20190730/202008_deep/20190730_(.*)_curated_to_R*"),
    
    # 20190812: mxMoM (separated folders per strain)
    myframes %>% ungroup %>% filter(date=="20200812",orientation=="b") %>% 
      group_by(gl_id,condition) %>% summarise(file_path=unique(file_path)) %>% 
      extract(file_path, "strain", "./data/Dany/20200812/20200812_bottom_(.*)_curated/20200812*"),
    
    # 20190812: mxMoM (separated folders per strain)
    myframes %>% ungroup %>% filter(date=="20200812",orientation=="t") %>% 
      group_by(gl_id,condition) %>% summarise(file_path=unique(file_path)) %>% 
      extract(file_path, "strain", "./data/Dany/20200812/20200812_top_(.*)_curated/20200812*"),
    
    myframes %>% ungroup %>% filter(date=="20200911",orientation=="b") %>% 
      group_by(gl_id,condition) %>% summarise(file_path=unique(file_path)) %>% 
      extract(file_path, "strain", "./data/Dany/20200911/20200911_bottom_(.*)_curated/20200911*"),
    
    myframes %>% ungroup %>% filter(date=="20200911",orientation=="t") %>% 
      group_by(gl_id,condition) %>% summarise(file_path=unique(file_path)) %>% 
      extract(file_path, "strain", "./data/Dany/20200911/20200911_top_(.*)_curated/20200911*")
  ))

myframes<- myframes %>% 
  separate(strain,c("vector","promoter"),remove=FALSE)

# CONVERT FLUO UNITS ####
#rmarkdown::render_site('./src/MoM_constitExpr_GFP_Estimation.Rmd')

myframes <- myframes %>%
  # append relevant conversion parameters
  mutate(autofluo_predict=NA, fp_per_dn=NA,
         # case of GFPmut2 (from zaslaver library)
         autofluo_predict = ifelse(strain!='asc662', 133.6, autofluo_predict),
         fp_per_dn = ifelse(strain!='asc662', 0.198, fp_per_dn),
         # case of asc662
         autofluo_predict = ifelse(strain=='asc662', 422.8, autofluo_predict),
         fp_per_dn = ifelse(strain=='asc662', 0.0361 * 4, fp_per_dn)) %>% 
  # convert to gfp units (after subtracting autofluorescence)
  mutate(fluogfp_amplitude = fluo_amplitude - autofluo_predict * length_um,
         gfp_nb = fluogfp_amplitude * fp_per_dn )
  
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
# Cell whose parents are known, that do not experience a switch, that divides, that are characterized by more than 4 data points

myframes_to_mycells <- myframes %>% 
  group_by(cell) %>% 
  filter(parent_id!=-1) %>% 
  filter(end_type=='div') %>% 
  filter(switch==FALSE) %>% 
  filter(n()>4) %>% 
  ungroup()

mycells <- myframes_to_mycells %>% 
  group_by(date,pos,orientation,gl,medium,vector,promoter,cell) %>% 
  partition(cluster=mycluster[1:nc]) %>%
  do((function(.df) {
    
    # .mod_ll_t <- lm( log(length_um)~time_sec, .df)  # use fastLm() for predict
    .mod_ll_t <- fastLmPure( cbind(1, .df$time_sec), log(.df$length_um) )
    .mod_lg_t <- fastLmPure( cbind(1, .df$time_sec), log(.df$gfp_nb) )
    .mod_l_t <- fastLmPure( cbind(1, .df$time_sec), .df$length_um )
    .mod_g_t <- fastLmPure( cbind(1, .df$time_sec), .df$gfp_nb )
    .mod_g_l <- fastLmPure( cbind(1, .df$length_um), .df$gfp_nb )
    
    .time_birth <- first(.df$time_sec)
    .time_div <- last(.df$time_sec)
    
    # .logl <- predict(.mod_ll_t, se.fit=TRUE)
    data.frame(npoints=.mod_ll_t$df.residual+1,
               time_birth=.time_birth, time_div=.time_div, div_time=.time_div-.time_birth,
               #cell_num_from_top=mean(.df$cell_num_in_lane),
               #cell_num_from_bottom=mean(.df$total_cell_in_lane-.df$cell_num_in_lane), 
               l_birth=first(.df$length_predict), l_div=last(.df$length_predict),
               # l_birth_raw=first(.df$length_raw),l_div_raw=last(.df$length_raw),
               # l_birth=exp(first(.logl$fit)), l_birth_se=exp(first(.logl$se.fit)), 
               # l_div=exp(last(.logl$fit)), l_div_se=exp(last(.logl$se.fit)), 
               logl_time_slope=.mod_ll_t$coefficients[2], logl_time_slopesd=.mod_ll_t$stderr[2],
               # logl_time_slope=.mod_ll_t$coefficients[2], logl_time_slopesd=summary(.mod_ll_t)$coefficients[2,2], 
               logl_time_r2=cor(.df$time_sec, log(.df$length_um))^2,
               logg_time_slope=.mod_lg_t$coefficients[2], logg_time_slopesd=.mod_lg_t$stderr[2], 
               logg_time_r2=cor(.df$time_sec, log(.df$gfp_nb))^2,
               l_time_slope=.mod_l_t$coefficients[2], l_time_slopesd=.mod_l_t$stderr[2], 
               l_time_r2=cor(.df$time_sec, .df$length_um)^2,
               g_birth=first(.df$gfp_nb), g_div=last(.df$gfp_nb), g_mean=mean(.df$gfp_nb),
               g_time_slope=.mod_g_t$coefficients[2], g_time_slopesd=.mod_g_t$stderr[2], 
               g_time_r2=cor(.df$time_sec, .df$gfp_nb)^2,
               g_l_slope=.mod_g_l$coefficients[2], g_l_slopesd=.mod_g_l$stderr[2], 
               g_l_r2=cor(.df$length_um, log(.df$gfp_nb))^2)
  })(.) ) %>% 
  collect() %>% 
  #arrange(condition, date, pos, gl, id) %>% 
  #mutate(gl_id=gsub('\\.[0-9]+$', '', cell)) %>% 
  ungroup()

mycells <- mycells %>% 
  # filter by r2 of exponential fit
  filter(logl_time_r2>0.95) %>% 
  # create new variables of interest
  mutate(dl=l_div - l_birth,
         alpha=log(l_div/l_birth) / div_time, # ok since l_xxx are fitted
         c_birth=g_birth/l_birth,
         c_div=g_div/l_div,
         dg=g_div - g_birth,
         dcdt=(g_div/l_div - g_birth/l_birth) / div_time,
         g=log(g_div/g_birth) / div_time, # would be cleaner from a fit but difficult to compare
         gamma=dg/dl,
         q=dg/dl*alpha)
  
myframes %>% distinct(cell) %>% nrow()
myframes_to_mycells %>% distinct(cell) %>% nrow()
mycells %>%  distinct(cell) %>% nrow()
mycells %>% group_by(promoter,medium) %>% count()
myframes_to_mycells <- myframes_to_mycells %>% 
  semi_join(mycells,by=c("cell"))



           