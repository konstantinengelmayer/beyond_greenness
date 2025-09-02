# Beyond Greenness

Code to reproduce the results of the study **“Beyond Greenness: Biophysical Index Groups and Their Ability to Explain Tree-Ring Growth Using a 23-Year Earth Observation Data Cube.”**

> All paths below are relative to the repository root.

## Reproduction

1) **Download SDC data** for the two tiles and place the files under `data/satellite_data/SDC/`.

- **Tiles**
  1. `kellerwald_lahntal` → files like `CSDC30_32UMB_YYYYDDD`
  2. `lindenberg_eifel_koenigsforst` → files like `CSDC30_31UGS_YYYYDDD`

- **Source**: <https://data-starcloud.pcl.ac.cn/aiforearth/#/>

- **Expected layout**
   ```
   data/
     satellite_data/
       SDC/
         kellerwald_lahntal/
           CSDC30_32UMB_YYYYDDD
           ...
         lindenberg_eifel_koenigsforst/
           CSDC30_31UGS_YYYYDDD
           ...
     vector_data/
       trees_all_plots.gpkg
       plot_metadata.xlsx
     dendro_data/
       MW3_Daten/
         Dendro-Data/
           d-XLS/
             ...
   ```

> `YYYYDDD` denotes year + day-of-year.

## Pipeline / Scripts

### 1) Rolling correlations (plot-wise)

**Run:**
```bash
Rscript 1_rolling_correlations_plotwise.r
```

**What it does:**

- Calls `extract_SDC_grouped_by_plot.r` to read SDC data from both `data/satellite_data/SDC/<tile>/` folders for each tree listed in `data/vector_data/trees_all_plots.gpkg`.
- Aggregates SDC by plot and writes:
  - `data/analysis_ready_data/SDC_extracted.csv`
- Computes vegetation index (VI) time series.
- Calls `extract_TRW.r` to:
  - Use the `tree_id` column in `data/vector_data/trees_all_plots.gpkg`
  - Extract tree-ring data per tree from `data/dendro_data/MW3_Daten/Dendro-Data/d-XLS/`
  - Detrend data and build plot-level TRI chronologies
- Matches TRI with SDC and runs the rolling correlation analysis.
- Saves the final correlation table to:
  - `data/analysis_ready_data/corr_table_vi_tri.rds`

> The dendro data used in this project can be obtained upon request.
### 2) Final figures

**Run:**
```bash
Rscript 2_final_figures.r
```

**What it does:**
- Computes all figures for the study’s paper and saves them to:
  - `figures/final_figures/`

## Additional data

The following helper scripts prepare auxiliary datasets:

- `a_sos_eos.r` — computes start and end of season (SOS/EOS) from SDC data.
- `a_topographics.r` — derives site topographic variables from SRTM30 data.
- `b_site_variables.r` — extracts climate variables for the sites from DWD climate data grids.

For precise details, see the header comments and inline documentation in each script.
