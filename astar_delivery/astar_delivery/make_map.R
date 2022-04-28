library(leaflet)
library(htmlwidgets)

df <- read.csv("results/spain_240949599_195977239.csv")
m <- leaflet(data = df) %>% addTiles() %>% addPolylines(lat = ~lat, lng = ~lon)
outname <- "results/spain_240949599_195977239.html"
title <- "A*_240949599_195977239.html"
saveWidget(m, file = outname, title = title)
