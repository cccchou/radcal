# resages Package

	resages is a package dedicated to the proper calculation of reservoir age offsets, particularly when dealing with pairs of radiocarbon (14C) and calendar ages. This package includes a set of radiocarbon calibration curves (IntCal13,  SHCal13，and IntCal20), a template folder, various examples, a manual, and three scripts written in python. To be precise, they are converted to python based on the Resage R language module(/ResAge_12_2015.zip)

## Installation

```sh
pip install resages
```

	**Note** Python needs to be **> version 3.10** because of the use of the `match` statement. If you don't want to upgrade Python, download the source code and change the `match` statement to an `if` statement.

## Usage

	Resages consists of three programs: `radical`, `rad2`, `colyear`. The corresponding function application scenarios can be found in the manual.This manual is the original author's source code manual for the R language

### Quick Use

```python
from resages.radical.radical import Radical
a = Radical("example_radcal_Southon2012", export_resage_pdf=True, threshold=1e-6)
a.radcal()

from resages.rad2.Rad2 import Rad2
a = Rad2(name="example_rad2")
a.rad2()

from resages.colyear.colyear import colyear
a = colyear(name="example_colyear")
a.colyear()
```


## Options
Each script has additional options:

`rad2`:

**calibrate**: Calibrates atmospheric 14C ages (default: FALSE).
**cc**: Selects the calibration curve (1 for IntCal20, 2 for Cal13, 3 for SHIntCal13).
**prob**: Confidence interval for calibrated ranges (default: 0.95).
**export_cal_pdf**: Exports the calibrated probability density functions (default: FALSE).

`colyear`:

**cc**: Selects the calibration curve (1 for IntCal20, 2 for Cal13, 3 for SHIntCal13).

`radcal`:

**mixture_pdf**: Computes the reservoir age offset probability density function mixture (default: FALSE).
**export_uncal_pdf**: Exports the uncalibrated probability density functions (default: FALSE).
**export_resage_pdf**:Exports the uncalibrated probability density functions (default: True, must set True because it export resage output.).
**cc**: Selects the calibration curve (1 for IntCal20, 2 for Cal13, 3 for SHIntCal13).

## Citation
This python program was adapted from the ##Resage R program## 
Please cite the following publications when using the ##ResAge package##:

- Soulet G, 2015. Methods and codes for reservoir-atmosphere 14C age offset calculations. Quaternary Geochronology 29:97-103, doi: 10.1016/j.quageo.2015.05.023
If using the F14R or δ14R metrics, cite: 
- Soulet G, Skinner LC, Beaupré SR, Galy V, 2016. A note on reporting of reservoir 14C disequilibria and age offsets. Radiocarbon, in press doi: 10.1017/RDC.2015.22

## References
- Beck, J.W., et al., 2001. Extremely large variations of atmospheric 14C concentration during the last glacial period. Science 292(5526): 2453-2458.
- Blaauw, M., 2010. Methods and code for ‘classical’ age-modelling of radiocarbon sequences. Quaternary Geochronology 5(5): 512-518.
- Bondevik, S., et al., 1999. Late Weichselian Marine 14C Reservoir Ages at the Western Coast of Norway. Quaternary Research 52(1): 104-114.
- Hall, B.L., et al., 2010. Constant Holocene Southern-Ocean 14C reservoir ages and ice-shelf flow rates. Earth and Planetary Science Letters 296(1): 115-123.
- Hoffmann, D.L., et al., 2010. Towards radiocarbon calibration beyond 28ka using speleothems from the Bahamas. Earth and Planetary Science Letters 289(1): 1-10.
- Hogg, A.G., et al., 2013. SHCal13 Southern Hemisphere calibration, 0–50,000 cal yr BP. Radiocarbon 55(4): 1889-1903.
Reimer, P.J., et al., 2009. IntCal09 and Marine09 Radiocarbon Age Calibration Curves, 0–50,000 years cal BP. Radiocarbon 51(4): 1111–1150.
- Reimer, P.J., et al., 2013. IntCal13 and Marine13 Radiocarbon Age Calibration Curves 0–50,000 Years cal BP. Radiocarbon 55(4): 1869-1887.
- Siani, G., et al., 2000. Radiocarbon reservoir ages in the Mediterranean Sea and Black Sea. Radiocarbon 42(2): 271-280.
Soulet, G., 2015. Methods and codes for reservoir-atmosphere 14C age offset calculations. Quaternary Geochronology 29: 97-103 doi: 10.1016/j.quageo.2015.05.023
- Southon, J., et al., 2012. A high-resolution record of atmospheric 14C based on Hulu Cave speleothem H82. Quaternary Science Reviews 33: 32-41.


