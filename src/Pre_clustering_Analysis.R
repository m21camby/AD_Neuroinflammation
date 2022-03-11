#!/bin RScript

library(Seurat)
library(ggplot2)
library(dplyr)
library(DT)
library(gridExtra)

sessionInfo()
library(magrittr)
library(devtools)
library(yaml)

# https://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
session_info() %>%
    write_yaml("./test.yaml")

