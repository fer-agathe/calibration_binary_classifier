This repository contains the R codes to replicate the results of the article titled _From Uncertainty to Precision: Enhancing Binary Classifier Performance through Calibration_.

The article is available on arXiv: <https://arxiv.org/abs/2402.07790>.

The replication codes can be found in the `replication` folder.

The `docs` folder contains an ebook with explanations of these codes. The HTML version of the ebook can be viewed at the following address: <https://fer-agathe.github.io/calibration_binary_classifier>.

# Replication Codes

The replication codes in the `replication` folder are divided in two parts, which can be run independently.

## Part 1: Synthetic Data and Calibration Metrics

In the first part, we provide codes to generate data from a know data generating process. We apply various calibration metrics and visualization techniques (file `sim-calibration.R`).

We then move to the presentation of the recalibration techniques and we apply them to our synthetic data (file `sim-recalibration.R`).

The functions contained in the `functions/simulations-synthetic.R` file are only used in our simulations. Therefore, they are not included in the `calibrationBinary` package.

R scripts required to reproduce the results:

- `sim-calibration.R`: calibration
- `sim-recalibration.R`: recalibration
- `functions/synthetic-data.R`: helper functions to generate data
- `functions/metrics.R`: helper functions to compute performance metrics and calibration metrics
- `functions/simulations-synthetic.R`: helper functions only used in this part.


## Part 2: Calibration of Random Forests with Real World Data

In this part, we use real-world data on default credit and investigate the calibration of random forests, as well as the relationship between calibration metrics and standard performance metrics.

R scripts:

- `rf-grid-search`: create datasets from data, train forests on a grid search
- `rf-simuls.R`: measure calibration
- `rf-simuls-result.R`: visualize the results.
- `functions/recalibration-methods.R`: recalibration function
- `functions/metrics.R`: helper functions to compute performance metrics and calibration metrics
- `functions/rf-functions`: other helper functions, specific to this second part.

# Reference of the paper

@misc{machado2024uncertainty,
   title={From Uncertainty to Precision: Enhancing Binary Classifier Performance through Calibration}, 
   author={Agathe Fernandes Machado and Arthur Charpentier and Emmanuel Flachaire and Ewen Gallic and François Hu},
   year={2024},
   eprint={2402.07790},
   archivePrefix={arXiv},
   primaryClass={cs.LG}
}
