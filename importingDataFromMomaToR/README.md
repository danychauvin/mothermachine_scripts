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


