import pandas as pd
import numpy as np
import matplotlib.pyplot
import seaborn
from pathlib import Path

DATA_DIR = Path("../data/raw/")
PATH_TO_CDG2019 = Path(f"{DATA_DIR}/cross_disorder/cdg2019/pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner.txt")

# TODO: Load first 1000 rows of the cdg2019 file
load_cdg2019_first_1000 = pd.read_csv(PATH_TO_CDG2019, delim_whitespace=True, nrows=1000)
