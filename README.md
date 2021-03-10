# ETCPy
**E**ffort-**T**o-**C**ompress in **Py**thon
 - [What is this](https://github.com/pranaysy/ETCPy#what-is-this)
   - [References](https://github.com/pranaysy/ETCPy#references)
 - [What can it do](https://github.com/pranaysy/ETCPy#what-can-it-do)
   - [Study Haemodynamics, Heart-Rate Variability and Cardiac Aging using ECG/EKG](https://github.com/pranaysy/ETCPy#study-haemodynamics-heart-rate-variability-and-cardiac-aging-using-ecgekg)
   - [Network Neuroscience, Psychophysics and Scientific Study of Consciousness](https://github.com/pranaysy/ETCPy#network-neuroscience-psychophysics-and-scientific-study-of-consciousness)
   - [Genome Complexity Analysis and Classification of Nucleotide Sequences](https://github.com/pranaysy/ETCPy#genome-complexity-analysis-and-classification-of-nucleotide-sequences)
   - [Audio Signal Processing and Denoising](https://github.com/pranaysy/ETCPy#audio-signal-processing-and-denoising)
 - [How to use it](https://github.com/pranaysy/ETCPy#how-to-use-it)
   - [Dependencies](https://github.com/pranaysy/ETCPy#dependencies)
   - [Installation](https://github.com/pranaysy/ETCPy#installation)
   - [Updating](https://github.com/pranaysy/ETCPy#updating)
   - [Usage](https://github.com/pranaysy/ETCPy#usage)
   - [Testing](https://github.com/pranaysy/ETCPy#testing)
   - [MATLAB Implementation](https://github.com/pranaysy/ETCPy#matlab-implementation)
 - [TODO](https://github.com/pranaysy/ETCPy#todo)
 - [License](https://github.com/pranaysy/ETCPy#license)

---

## What is this
A Python implementation of the compression-complexity measure called Effort-To-Compress or ETC. ETC captures the compressibility and complexity of discrete symbolic sequences using lossless compression. It has been shown to robustly estimate complexity, comparing favorably for short and noisy time series in comparison with entropy and Lempel-Ziv complexity.

Using ETC, causal information flow between multiple discrete symbolic sequences can be assessed and recently, such a use has been presented, rigorously proven and demonstrated to be an effective model-free measure of causality. Introduced as Compression-Complexity Causality or CCC, this measure is robust to numerous data contaminants, noise sources and pre-processing artifacts. On comparison with Granger Causality and Transfer Entropy, CCC compares favorably and outperforms them on synthetic as well as real world causal interactions. An implementation of CCC is included in this repository.

While any lossless compressor may be used with ETC and subsequently with CCC, a grammar-based lossless compression algorithm called Non-Sequential Recursive Pair Substitution or NSRPS is used presently. NSRPS has been rigorously studied and shown to be an effective tool for data compression and entropy estimation. This repository also contains a fast Cython implementation of NSRPS for use with ETC and CCC.

#### References
 - Benedetto, Dario, Emanuele Caglioti, and Davide Gabrielli. “Non-Sequential Recursive Pair Substitution: Some Rigorous Results.” Journal of Statistical Mechanics: Theory and Experiment 2006, no. 09 (September 25, 2006): P09011–P09011. https://doi.org/10.1088/1742-5468/2006/09/P09011.
 - Balasubramanian, Karthi, Gayathri R. Prabhu, Lakshmipriya V. K. , Maneesha Krishnan, Praveena R. , and Nithin Nagaraj. “Classification of Periodic, Chaotic and Random Sequences Using NSRPS Complexity Measure.” ArXiv:1205.4886 [Nlin], May 22, 2012. http://arxiv.org/abs/1205.4886.
 - Nagaraj, Nithin, Karthi Balasubramanian, and Sutirth Dey. “A New Complexity Measure for Time Series Analysis and Classification.” The European Physical Journal Special Topics 222, no. 3–4 (July 2013): 847–60. https://doi.org/10.1140/epjst/e2013-01888-9.
 - Nagaraj, Nithin, and Karthi Balasubramanian. “Dynamical Complexity of Short and Noisy Time Series: Compression-Complexity vs. Shannon Entropy.” The European Physical Journal Special Topics 226, no. 10 (July 2017): 2191–2204. https://doi.org/10.1140/epjst/e2016-60397-x.
 - Kathpalia, Aditi, and Nithin Nagaraj. “Data-Based Intervention Approach for Complexity-Causality Measure.” PeerJ Computer Science 5 (May 27, 2019): e196. https://doi.org/10.7717/peerj-cs.196.


## What can it do
#### Study Haemodynamics, Heart-Rate Variability and Cardiac Aging using ECG/EKG
   - Balasubramanian, Karthi, Nithin Nagaraj, and Sandipan Pati. “Chaos or Randomness? Effect of Vagus Nerve Stimulation During Sleep on Heart-Rate Variability.” IETE Journal of Research, June 30, 2020, 1–7. https://doi.org/10.1080/03772063.2020.1780165.
   - Srilakshmi, P, Karthi Balasubramanian, Nithin Nagaraj, and Sandipan Pati. “Multiscale Analysis of Heart Rate Variability Using Subsymmetry and Effort-to-Compress Complexity Measures.” In 2018 15th IEEE India Council International Conference (INDICON), 1–5. Coimbatore, India: IEEE, 2018. https://doi.org/10.1109/INDICON45594.2018.8986972.
   - Thanaj, Marjola, Andrew J. Chipperfield, and Geraldine F. Clough. “Analysis of Microvascular Blood Flow and Oxygenation: Discrimination between Two Haemodynamic Steady States Using Nonlinear Measures and Multiscale Analysis.” Computers in Biology and Medicine 102 (November 2018): 157–67. https://doi.org/10.1016/j.compbiomed.2018.09.026.
   - Balasubramanian, Karthi, K Harikumar, Nithin Nagaraj, and Sandipan Pati. “Vagus Nerve Stimulation Modulates Complexity of Heart Rate Variability Differently during Sleep and Wakefulness.” Annals of Indian Academy of Neurology 20, no. 4 (2017): 403. https://doi.org/10.4103/aian.AIAN_148_17.
   - Balasubramanian, Karthi, and Nithin Nagaraj. “Aging and Cardiovascular Complexity: Effect of the Length of RR Tachograms.” PeerJ 4 (2016): e2755. https://doi.org/10.7717/peerj.2755.

#### Network Neuroscience, Psychophysics and Scientific Study of Consciousness
   - Ashley J. Funkhouser. "The Role of Action in Affordance Perception Using Virtual Reality" 2020. Honors College Thesis with Dr. Alen Hajnal, Department of Psychology, The University of Southwestern Mississipi. https://aquila.usm.edu/honors_theses/714/
   - Agarwal, Nikita, Aditi Kathpalia, and Nithin Nagaraj. “Distinguishing Different Levels of Consciousness Using a Novel Network Causal Activity Measure.” In 2019 Global Conference for Advancement in Technology (GCAT), 1–5. BANGALURU, India: IEEE, 2019. https://doi.org/10.1109/GCAT47503.2019.8978424.
   - Virmani, Mohit, and Nithin Nagaraj. “A Novel Perturbation Based Compression Complexity Measure for Networks.” Heliyon 5, no. 2 (February 2019): e01181. https://doi.org/10.1016/j.heliyon.2019.e01181.
   - Kondo, Fumika. “Can Alterations in the Temporal Structure of Spontaneous Brain Activity Serve as a Disease-Specific Biomarker for Schizophrenia? A Multi Cohort FMRI Study,” 2017. https://doi.org/10.20381/RUOR-20801.
   - Kimiskidis, Vasilios K., Christos Koutlis, Alkiviadis Tsimpiris, Reetta Kälviäinen, Philippe Ryvlin, and Dimitris Kugiumtzis. “Transcranial Magnetic Stimulation Combined with EEG Reveals Covert States of Elevated Excitability in the Human Epileptic Brain.” International Journal of Neural Systems 25, no. 05 (August 2015): 1550018. https://doi.org/10.1142/S0129065715500185.

#### Genome Complexity Analysis and Classification of Nucleotide Sequences
   - Balasubramanian, Karthi, and Nithin Nagaraj. “Automatic Identification of SARS Coronavirus Using Compression-Complexity Measures.” Preprint. Bioinformatics, March 27, 2020. https://doi.org/10.1101/2020.03.24.006007.

#### Audio Signal Processing and Denoising
   - Kiefer, Chris, Overholt, Dan and Eldridge, Alice (2020) Shaping the behaviour of feedback instruments with complexity-controlled gain dynamics. New Interfaces for Musical Expression, Birmingham, UK, 21-25 July 2020. Published in: Proceedings of the International Conference on New Interfaces for Musical Expression. 343-348. NIME, Birmingham, UK. ISSN 2220-4806. https://sro.sussex.ac.uk/id/eprint/91009/
   - Li, Guohui, Qianru Guan, and Hong Yang. “Noise Reduction Method of Underwater Acoustic Signals Based on CEEMDAN, Effort-To-Compress Complexity, Refined Composite Multiscale Dispersion Entropy and Wavelet Threshold Denoising.” Entropy 21, no. 1 (December 24, 2018): 11. https://doi.org/10.3390/e21010011.


## How to use it
The simplest way right now is to use `pip` to clone this repository and install locally inside a `conda` or a `virtualenv` environment. This way several functions implemented in Cython will be automatically compiled natively on the host system. Instructions below.

While the repository is called `ETCPy`, the package namespsace available for use is `ETC`. All functionality is available through the `ETC` namespace.

For running tests (strongly recommended), additional packages need to be installed.

### Operating System Support
 - GNU/Linux-based distributions (tested on Ubuntu 16.04, 18.04, 20.04)
 - **Currently does not work out of the box on Windows.** Cython and C/C++ build toolchain need to be setup properly for compilation on Windows to work. It may work with some gymnastics using MinGW + Visual Studio Build Tools, **currently untested.** Although does work on WSL!

### Dependencies
For core functionality:
 - `numpy`
 - `pandas`
 - `joblib`
 - `cython`
   - Note: Cython needs a working C/C++ compiler such as GCC/Clang and associated build-utils/toolchain. While it should work out of the box on any modern Linux distribution, ensure a proper installation as instructed in the [official documentation.](https://cython.readthedocs.io/en/latest/src/quickstart/install.html).

For tests:
 - `pytest`
 - `hypothesis`

### Installation
Skip the first step if an environment is already available:
1. Create a fresh `conda` or `pip`/`virtualenv`-based environment with `numpy` and `cython` packages. Choose an appropriate name instead of `myenv`.
    ```bash
    $ conda create -n myenv python numpy pandas joblib cython
    ```
2. Activate environment using `conda activate myenv` or virtualenv equivalent.

   If `git` is not installed, then:
     - either install it at a system level directly from the [official website](https://git-scm.com/download) or via prefereed package manager
     - or install it within the newly created conda environment using `conda install git`

3. Use `pip`* to install directly from GitHub using the `git` VCS backend
    ```bash
    $ python -m pip install git+https://github.com/pranaysy/ETCPy.git
    ```
4. Done! Open a Python shell, execute `import ETC` and proceed to the [demo](./demo.py)

---
*mixing `pip` and `conda` is not a generally advised but can be used based on [certain recommendations](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#pip-in-env)

### Updating
Use the `-U` flag with pip for updating to the most current version available from this repository:
```
$ python -m pip install -U git+https://github.com/pranaysy/ETCPy.git
```
This will rebuild the compiled Cython functions as well.

### Usage
Please check out [`demo.py`](./demo.py) to see ETC in action. [Functions for dealing with NumPy arrays](https://github.com/pranaysy/ETCPy/blob/master/demo.py#L121) are also available. In addition to the core functionality of ETC, a [brief demo of Compression-Complexity Causality (CCC)](https://github.com/pranaysy/ETCPy/blob/master/demo.py#L158) is also included for uncoupled as well as coupled first-order auto-regressive processes.

The implementations of ETC as well as CCC include multicore parallelization (using [`joblib`](https://joblib.readthedocs.io/en/latest/index.html)) and can benefit from more available CPU cores for multiple sequences.

### Testing
Most of the tests are property-based or behavior-based, and are implemented using the awesome [`hypothesis` framework](https://hypothesis.readthedocs.io/en/latest/).
Make sure dependencies are satisfied within the working environment:
```bash
$ python -m pip install -U pytest hypothesis
```
Grab a copy of this repository using git and enter the local directory:
```bash
$ git clone https://github.com/pranaysy/ETCPy.git
$ cd ETCPy
```
Run tests:
```bash
$ pytest ETC/
```

### MATLAB Implementation
 - The original ETC implementation in MATLAB can be found here: https://sites.google.com/site/nithinnagaraj2/journal/etc


## TODO
 - Hyperparameter optimization for CCC
 - Add performance metrics
 - Automated tests with `tox`
 - Better packaging: `pip` vs `conda`
 - Visualizations
 - Improve test coverage
 - Documentation using Sphinx/MkDocs
 - Windows support

## License
Copyright 2021 Pranay S. Yadav and Nithin Nagaraj

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
