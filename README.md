sclero
======
**An R package to measure growth patterns and align sampling spots in photographs**

This is the developmental version of the *sclero* package providing functions to measure growth patterns and align sampling spots in chronologically deposited materials. The package is intended for the fields of sclerochronology, dendrochronology and geology. A stable version of the package can be found from [CRANsclero][CRAN].

Installation of the GitHub version
-------
Using the *devtools* package:
```{r}
library(devtools)
install_github("MikkoVihtakari/sclero", dependencies = TRUE)
```

User Manual
-------
A pdf version of the user manual is available under [`/inst/doc`][doc] and a source version under [`/vignettes`][vignettes] as an .Rnw file.

Contributions and contact information
-------
Any contributions to the package are more than welcome. Please contact the package creator Mikko Vihtakari (<mikko.vihtakari@gmail.com>) to discuss your ideas on improving the package.

Dependencies
--------
The *sclero* package depends on:
- [RImageJROI][RImageJROI] package. Used to import ImageJ ROI objects to R.
- [spatstat][spatstat] package. Used for geometric calculations.
- [plyr][plyr] package. Used for quicker and easier list calculations. This dependency could potentially be removed, but the removal would require rewriting many functions included in the *sclero* package

[CRANsclero]: 
[RImageJROI]: https://github.com/davidcsterratt/RImageJROI
[spatstat]: http://cran.r-project.org/web/packages/spatstat/index.html
[plyr]: https://github.com/hadley/plyr
[doc]: https://github.com/MikkoVihtakari/sclero/tree/master/inst/doc
[vignettes]: https://github.com/MikkoVihtakari/sclero/tree/master/vignettes
