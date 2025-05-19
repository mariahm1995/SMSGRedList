
**Tools for spatial reassessment of small mammal species for the IUCN Red List**, developed for use by the IUCN SSC Small Mammal Specialist Group.\
\
---\
\
## What is SMSGRedList?\
\
`SMSGRedList` is an R package that automates and standardizes steps in the spatial component of Red List reassessments. It supports:\
\
- Downloading and cleaning GBIF and literature occurrence data\
- Generating new range polygons (MCP, kernel, alpha hull)\
- Applying land, elevation, and manual corrections\
- Calculating EOO and AOO\
- Producing expert-review maps (leaflet format)\
- Exporting IUCN-ready shapefiles\
\
---\
\
## Requirements\
\
Before using this package, you **must download a folder of base spatial data** hosted on **Figshare**. This folder contains:\
\
- Country boundaries\
- Elevation rasters\
- IUCN polygon templates\
- 2\'d72 km AOO grids\
\
\uc0\u55357 \u56599  **Figshare download link**:  \
https://figshare.com/account/projects/247241/articles/28915925?file=54137138\
\
Once downloaded, place the folder in your project directory as "data"}
