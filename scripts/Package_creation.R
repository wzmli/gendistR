# Package creation notes

# from https://r-pkgs.org/whole-game.html

library(devtools)
use_mit_license()
load_all()
document()
check()
install()
library(gendistR)
?gendist_avg
?gendist_parse
build_readme()
usethis::use_package("dplyr") # add dependencies to imports in description
