# Sun Induced Fluroscence Downscaling processor using openEO and Copernicus Data Space Ecosystem Infrastructure

## Description


TBD


## Tutorial

You can run the processor following the [tutorial](https://dpabon.github.io/SIF_downscaling_CDSE/)

- First you need to install conda-forge:

- clone this github repository:

```git clone https://github.com/dpabon/SIF_downscaling_CDSE```

- Then re-create the conda environment:

```cd SIF_downscaling_CDSE```
```conda env create -f environment.yml```

Now you have everything setup to run the tutorial.


- ```openEO_sif_downscaling.py``` contains the active development of the SIF downscaling workflow using openEO.
- ```udf.py``` contains the User Defined Function need it for openEO.
- ```environment.yml``` contains the conda environment with all the packages need it to reproduce the analysis.
- ```data``` contains a COG file and the corresponding geojson STAC.


## FAQ

- Q: I have problems to run the tutorial?
    A: Please open an [issue](https://github.com/dpabon/SIF_downscaling_CDSE/issues/new) 

- Q: Can I select a new area to apply the procesor?
 A: Sure, just change the values in the area of interest cell

- Q: How can I contribute?
 A: clone this repository and create pull requests.


## Acknowledgement


This project has received funding from the [Open-Earth-Monitor Cyberinfrastructure](https://earthmonitor.org/) project that is part of European Union's Horizon Europe research and innovation programme under grant [101059548](https://cordis.europa.eu/project/id/101059548).