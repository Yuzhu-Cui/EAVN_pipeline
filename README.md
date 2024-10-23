# EAVN_pipeline
This pipeline is primarily designed for the data reduction of the East Asia VLBI Network ([EAVN](https://radio.kasi.re.kr/eavn/main.php)) but is versatile enough to be used for general VLBI data processing. The pipeline automates the entire data reduction process, integrating command-line tools for tasks such as data ingestion, editing, calibration, and imaging into a streamlined workflow. The process is divided into two key stages: calibration and imaging. In the calibration stage, the pipeline handles the preparation and refinement of the data, ensuring accurate and high-quality input for the imaging stage, where the final science-ready images are produced. This automation reduces human error and increases processing efficiency for large datasets.


## Dependencies
- [PGPLOT](https://sites.astro.caltech.edu/~tjp/pgplot/) (version=5.2)
- [AIPS](http://www.aips.nrao.edu/index.shtml) (version=31DEC24)
- [Parseltongue](https://www.jive.eu/jivewiki/doku.php?id=parseltongue:parseltongue) (version=3.0)
- [Difmap](https://science.nrao.edu/facilities/vlba/docs/manuals/oss2013a/post-processing-software/difmap) (version=2.5q)
- [Anaconda](https://www.anaconda.com/) (version=2021.05)  
Detailed installation steps can be found in: [install.md](https://github.com/lao19881213/EAVN_pipeline/blob/main/install.md).

## Usage
### Set up software environment  
- Add all dependency environments to a file (for example, bashrc), and then: 
```
source bashrc
```
### Single data file processing 
- Calibration   
1. Edit the input file (template.inp) and set the required parameters  
2. Run the data calibration pipeline:  
```
ParselTongue EAVN.py template.inp
```
- Imaging  
1. Edit the input file (imaging.inp) and set the required parameters  
2. Run the imaging pipeline:
```
python3 difmap_imaging.py imaging.inp
``` 
### Multiple data files batch processing  
```
ParselTongue batch_processing.py  
```

## Contributors
Baoqiang Lao, Yuzhu Cui

## References
- [Cui, Y., Hada, K., on behalf of EAVN Science Working Group, Data reduction memo: Amplitude calibration guideline of TMRT and NSRT (2020).](https://radio.kasi.re.kr/eavn/pdf/Amplitude_calibration_guideline_of_TMRT_and_NSRT.pdf)
- [Cui, Y., Hada, K. et al. East Asian VLBI Network observations of active galactic nuclei jets: imaging with KaVA+Tianma+Nanshan. Res. Astron. Astrophys. 21, 205 (2021).](https://www.raa-journal.org/issues/all/2021/v21n8/202203/t20220323_21975.html)
- [Cui, Y., Hada, K., Kawashima, T. et al. Precessing jet nozzle connecting to a spinning black hole in M87. Nature 621, 711–715 (2023).](https://www.nature.com/articles/s41586-023-06479-6)

## Acknowledgements
We thank Kazuhiro Hada for his valuable guidance in the data calibration process. The pipeline code was developed based on the [EVN pipeline](https://www.jive.eu/jivewiki/doku.php?id=parseltongue:grimoire). Originally crafted by [Cormac Reynolds](reynolds@jive.nl), the EVN pipeline is now maintained and further developed by [Stephen Bourke](bourke@jive.nl). For those seeking a concise guide on executing the pipeline, a brief manual is accessible at [pypeline_public.pdf](https://www.jive.eu/jivewiki/lib/exe/fetch.php?media=parseltongue:pypeline_public.pdf). This project is supported by the China Postdoctoral Science Foundation (No. 2024T170845) and the Natural Science Foundation of China (grant 12303021).
