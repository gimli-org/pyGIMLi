Reproducibility tips (for beginners)
====================================

• Use a pinned environment (conda):
  conda create -n pg -c gimli -c conda-forge pygimli
  conda activate pg

• Set random seeds in examples to stabilize results:
  import numpy as np; np.random.seed(0)

• Note OS differences (Windows vs. Linux/Mac) when paths or installers differ.

• Link to install & contributing guides for more detail.

References: Installation, Contributing.