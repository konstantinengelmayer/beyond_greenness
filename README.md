# Beyond Greenness
This repo contains the code to reproduce the findings of the study "Beyond Greenness: Biophysical Index Groups and Their Ability to Explain Tree-Ring Growth Using a 23-Year Earth Observation Data Cube"

## Reproduction
To reproduce the results of the study first the SDC data of the two tiles `kellerwald_lahntal` and `lindenberg_eifel_koenigsforst` needs to be downloaded and saved in the correct folders in the `data/satellite_data/SDC` folder. The data can be downloaded here: <https://data-starcloud.pcl.ac.cn/aiforearth/#/>

#### Tile names:
kellerwald_lahntal: `CSDC30_32UMB_yyyyddd`
lindenberg_eifel_koenigsforst: `CSDC30_31UGS_yyyyddd`

## Scripts
1. First the script `1_rolling_correlations_plotwise.r` needs to be run. The script runs the function `extract_SDC_grouped_by_plot.r` and extracts the SDC data found in both `data/satellite_data/SDC` folders for each tree of the file `data/vector_data/trees_all_plots.gpkg`. The extracted file is saved as `data/analysis_ready_data/corr_table_vi_tri.rds`
