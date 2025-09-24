# Data package: GIS files of nocturnal bird migration over the Netherlands

This data package accompanies the manuscript:
Large-scale mapping of nocturnal bird migration to accelerate a nature-inclusive energy transition

## Authors (and affiliations)

Bart Hoekstra (1), Bart Kranstauber (1), Maja Bradarić (1), Johannes De Groeve (1), Stacy Shinneman (1), Berend C. Wijers (1), Hidde Leijnse (3), Hans van Gasteren (1,4), Adriaan M. Dokter (1,2), Emiel van Loon (1), Judy Shamoun-Baranes (1)

1. Institute for Biodiversity and Ecosystem Dynamics, University of Amsterdam, P.O. Box 94240, 1090 GE Amsterdam, The Netherlands<br />
2. Cornell Lab of Ornithology, Cornell University, 159 Sapsucker Woods Rd, Ithaca, NY 14850, United States of America<br />
3. R&D Observations and Data Technology, Royal Netherlands Meteorological Institute, De Bilt, The Netherlands<br />
4. Royal Netherlands Air Force, P.O Box 8762, 4820 BB Breda, The Netherlands<br />

# Data

The data package contains composite maps derived from Dutch weather radar data, summarizing peak nights of nocturnal bird migration during spring and autumn seasons.

The maps are provided as raster files (.tif) at 500 × 500 m resolution. An accompanying .geojson file delineates areas of uncertainty (e.g. from ground clutter and beam blockage). This file should always be used together with the rasters when interpreting the data.

For details on the methodology and guidance on interpretation, please refer to the associated paper.

If you use these data, please cite the paper (see repository/Zenodo entry for full citation).

## File overview

Filename	                        Description
NL_autumn_seasonal_passage.tif    Seasonal passage map for autumn (Figure 2 in the paper)
NL_autumn_bird_density.tif        Final map of bird densities in autumn
NL_autumn_VIR_raw.tif             Raw (unsmoothed) vertically integrated reflectivity map
NL_spring_seasonal_passage.tif    Seasonal passage map for spring (Figure 2 in the paper)
NL_spring_bird_density.tif        Final map of bird densities in spring
NL_spring_VIR_raw.tif             Raw (unsmoothed) vertically integrated reflectivity map
NL_uncertainty_areas.geojson      Areas affected by ground clutter or beam blockage

## Links
- GitHub repository with code and workflow: https://github.com/barthoekstra/nl-birdmap
- Archived version on Zenodo (with DOI): [link once available]
