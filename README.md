# nl-birdmap
This repository accompanies the manuscript:

_Large-scale mapping of nocturnal bird migration to accelerate a nature-inclusive energy transition_

It contains the code and data necessary to reproduce the key analyses and figures.


## Authors

Bart Hoekstra<sup>1</sup>, Bart Kranstauber<sup>1</sup>, Maja Bradarić<sup>1</sup>, Johannes De Groeve<sup>1</sup>, Stacy Shinneman<sup>1</sup>, Berend C. Wijers<sup>1</sup>, Hidde Leijnse<sup>3</sup>, Hans van Gasteren<sup>1,4</sup>, Adriaan M. Dokter<sup>1,2</sup>, Emiel van Loon<sup>1</sup>, Judy Shamoun-Baranes<sup>1</sup>

<sup>1</sup> Institute for Biodiversity and Ecosystem Dynamics, University of Amsterdam, P.O. Box 94240, 1090 GE Amsterdam, The Netherlands<br />
<sup>2</sup> Cornell Lab of Ornithology, Cornell University, 159 Sapsucker Woods Rd, Ithaca, NY 14850, United States of America<br />
<sup>3</sup> R&D Observations and Data Technology, Royal Netherlands Meteorological Institute, De Bilt, The Netherlands<br />
<sup>4</sup> Royal Netherlands Air Force, P.O Box 8762, 4820 BB Breda, The Netherlands<br />

## Structure of this repository
This repository is organized around three key scripts:
1. `radar-profiles.Rmd`: Derives flight altitudes from vertical profile radar data. Produces Figure 1 of the manuscript.
2. `birdmap.Rmd`: Contains the processing from identified peak nights to the peak night distributions of migratory birds across the Netherlands, and the subsequent analysis of wind power-bird migration tradeoffs. Produces figure 2 and 3 of the manuscript.
3. `rbc_preprocessing.R`: Performs computationally intensive pre-processing of radar data to correct for range bias. Best run independently before executing `birdmap.Rmd`.


## Data
This repository contains all data needed as initial input data, or GIS data to recreate the figures. Unfortunately, given the large data volume of weather radar data, it is impractical to share this here, but the data is freely available via the [KNMI open data portal](https://dataplatform.knmi.nl/) and can be accessed via the new [{getRad} R package](https://aloftdata.github.io/getRad/). The analysis notebooks show how files should be structured for the analysis and the folder structure has been retained by adding empty `.gitkeep` files to this repository.

A short summary table with of input files, their contents and sources is listed below. All files are located in `data` and the respective subfolders. A dash indicated in the source column indicates the data is created for the purpose of this study.

| Filename | Subfolder | Source | Contents | URL | Comment |
|----------|--------|--------|----------|-----|---------|
|`vp_mtr_ordered_seasonal_day_`<br />`night_migration.csv`|        |-        |Ranked nightly migration activity across radars derived from vertical profiles of birds     |     |Main input file for the study, also contains information on daytime migration|
|`dhl.csv`          |`beam_blockage`        |-       |Beam blockage for Den Helder radar          |     |Generated with `beam-blockage/beamblockage.ipynb`          |
|`hrw.csv` | `beam_blockage` |-|Beam blockage for Herwijnen radar | |Generated with `beam-blockage/beamblockage.ipynb`|
|`NL_Beamblockage.geojson`| `gis` |-|Beam blockage for Dutch radars converted to .geojson|Derived from above|
|`manual_clutter.geojson` | `gis` |-|Manually drawn areas affected by ground clutter | |Only if not resolved by beam blockage procedure|
|`NL_Clutter-500m-Buffer.geojson`| `gis` |-|.geojson of identified clutter areas buffered by 500m| |
|`Netherlands-5km-Buffer.geojson`| `gis` |-|.geojson with an outline of the Netherlands buffered by 500m| |
|`Netherlands.geojson`| `gis` |-|.geojson with an outline of the Netherlands| |
|`Noord-Holland-10km-Buffer.geojson`| `gis` |-|.geojson with an outline of the province of North-Holland bufffered by 10km| |
|`Noord-Holland.geojson`| `gis` |-|.geojson with an outline of the province of North-Holland| |
|`Radars-100km-Buffer.geojson`| `gis` |-|.geojson with the area enclosed within 100km from both Dutch weather radars| |
|`Radars.geojson`| `gis` |-|.geojson with the locations from both Dutch weather radars| |
|`gemeenten.geojson`| `gis` | CBS |.geojson with the Dutch municipalities| [URL](https://www.cbs.nl/nl-nl/dossier/nederland-regionaal/geografische-data/cbs-gebiedsindelingen)| Not in final version of manuscript|
|`rivm_20250101_windturbines_`<br />`ashoogte.shp` | `gis` | RIVM|Shapefile with turbine locations and dimension information|[URL](https://www.nationaalgeoregister.nl/geonetwork/srv/api/records/23d0d402-a6d9-47c5-a6f3-d7f7fb35cb79?language=all) | Not in final version of manuscript|
|`altitude_profiles_weatherradar.rds`| `profiles` |-|Aggregated altitude profiles from both Dutch weather radars| |
|`altitude_profiles_weatherradar.sql`| `profiles` |-|SQL Query to aggregate altitude profiles from both Dutch weather radars| |
|`RES_zoekgebieden_*`|`gis/RES/`|Province of North-Holland|Shapefiles with information on the search areas/candidate sites for renewable energy developments|[URL](https://apps.vertigisstudio.eu/web/?app=194abc647f794375873dcd563932dd8e)|
|`vp_mtr_cum50_cl.csv`||-|.csv export of vp_mtr_cum50_cl object containing scan screening outcomes||
