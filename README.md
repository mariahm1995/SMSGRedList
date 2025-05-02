{\rtf1\ansi\ansicpg1252\cocoartf2761
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # SMSGRedList <img src="https://raw.githubusercontent.com/rstudio/hex-stickers/master/PNG/Rlogo.png" height="25" />\
\
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
Once downloaded, place the folder in your project directory as a \'93data\'94}