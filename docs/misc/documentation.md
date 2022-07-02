# Generate documentation

## How to generate documentation

Technologies used: Doxygen + Spnynx + Breathe

to be written ...

## Create pdf documentation from Rd files

* Install volesti library.
* In `R` mode (or in Rstudio) Run
```r
pack = "volesti"
path = find.package(pack)
system(paste(shQuote(file.path(R.home("bin"), "R")),
    "CMD", "Rd2pdf", shQuote(path)))
```
* The pdf will be created and saved in `R-proj` folder.
* We give such a documentation in `/R-proj/doc` folder.
