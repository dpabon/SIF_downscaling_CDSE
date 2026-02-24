# Sun Induced Fluroscence Downscaling processor using openEO and Copernicus Data Space Ecosystem Infrastructure

## Description


Sun-induced chlorophyll fluorescence (SIF) is a re-emitted signal directly originating from the photosynthetic system, representing plant fluorescence observed from space. It is therefore considered a direct measure of plant physiological status and responses to environmental changes, while also serving as a potential indicator for estimating global gross primary productivity (GPP). However, the revisit cycles and operational altitudes of existing SIF observation platforms limit their spatio-temporal resolution. This results in decoupling between the linear relationship between global-scale SIF and GPP at finer temporal and spatial scales.
The Sentinel-5P TROPOMI sensor, with its daily nadir spatial resolution of 3.5 × 5.5-7 km, enables more detailed observation of SIF/GPP in terrestrial ecosystems. The following example illustrates the downscaling process (i.e., increasing the spatial resolution) for SIF based on the TROPOMI SIF.


## Tutorial

You can run the processor following the [tutorial](https://dpabon.github.io/SIF_downscaling_CDSE/) (We recommend [Positron](https://positron.posit.co/) as it included all the batteries needed) 

- First you need to install pixi:

https://pixi.sh/latest/

- In positron set pixi tool path:

To check where pixi was installed you can run in a terminal

Linux and Mac:
```whereis pixi```

Windows:
```where pixi```

and copy paste the path into Positron settings "Python: Pixi Tool Path".

- You need to fork this github repository on github (This step is necessary as currently CDSE is not able to perform upsample operations. Then a spatial upsampling needs to run locally, upload to github and load the results in the CDSE again):

![fork figure](fork.png)

Then clone the repo into your local machine:

```git clone https://github.com/your_user_name/SIF_downscaling_CDSE```

- Then install the dependencies using pixi:

```cd SIF_downscaling_CDSE```

```pixi install```

Now you have everything setup to run the [tutorial](https://dpabon.github.io/SIF_downscaling_CDSE/), Don't forget to select the pixi python interpreter before running.


- ```openEO_sif_downscaling.py``` contains the active development of the SIF downscaling workflow using openEO.
- ```udf.py``` contains the User Defined Function need it for openEO.
- ```environment.yml``` contains the conda environment with all the packages need it to reproduce the analysis.
- ```data``` contains a COG file and the corresponding geojson STAC.


## FAQ

- Q: Problems to run the tutorial?

A: Please open an [issue](https://github.com/dpabon/SIF_downscaling_CDSE/issues/new) 

- Q: Can I select a new area to apply the procesor?
 
A: Sure, just change the values in the area of interest cell

- Q: Why pixi if everyone is using coda-forge?

A: Mainly because of convenience. pixi allows to easily create workspace for multiple platforms (the ones for this project include "linux-64", "win-64", "osx-64", "osx-arm64").

- Q: How can I contribute?

A: clone this repository and create pull requests.


## Acknowledgement


This project has received funding from the [Open-Earth-Monitor Cyberinfrastructure](https://earthmonitor.org/) project that is part of European Union's Horizon Europe research and innovation programme under grant [101059548](https://cordis.europa.eu/project/id/101059548).