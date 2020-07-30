# ETCPy

## What is this
A Python implementation of the compression-complexity measure called Effort-To-Compress or ETC. ETC captures the compressibility and complexity of discrete symbolic sequences using lossless compression. It has been shown robustly estimate complexity, comparing favorably for relatively short and noisy time series relative to entropy and Lempel-Ziv complexity.

Using ETC, causal information flow between multiple discrete symbolic sequences can be assessed and recently, such a use has been presented, rigorously proven and demonstrated to be an effective model-free measure of causality. Introduced as Compression-Complexity Causality or CCC, this measure is robust to numerous data contaminants, noise sources and pre-processing artifacts. On comparison with Granger Causality and Transfer Entropy, CCC compares favorably and outperforms them on synthetic as well as real world benchmarks. An implementation of CCC is included in this repository.

While any lossless compressor may be used with ETC and subsequently with CCC, a grammar-based lossless compression algorithm called Non-Sequential Recursive Pair Substitution or NSRPS is used presently. NSRPS has been rigorously studied and shown to be an effective tool for data compression and entropy estimation. This repository also contains a fast Cython implementation of NSRPS for use with ETC and CCC.

#### References
 - Balasubramanian, Karthi, Gayathri R. Prabhu, Lakshmipriya V. K. , Maneesha Krishnan, Praveena R. , and Nithin Nagaraj. “Classification of Periodic, Chaotic and Random Sequences Using NSRPS Complexity Measure.” ArXiv:1205.4886 [Nlin], May 22, 2012. http://arxiv.org/abs/1205.4886.
 - Benedetto, Dario, Emanuele Caglioti, and Davide Gabrielli. “Non-Sequential Recursive Pair Substitution: Some Rigorous Results.” Journal of Statistical Mechanics: Theory and Experiment 2006, no. 09 (September 25, 2006): P09011–P09011. https://doi.org/10.1088/1742-5468/2006/09/P09011.
 - Kathpalia, Aditi, and Nithin Nagaraj. “Data-Based Intervention Approach for Complexity-Causality Measure.” PeerJ Computer Science 5 (May 27, 2019): e196. https://doi.org/10.7717/peerj-cs.196.
 - Nagaraj, Nithin, Karthi Balasubramanian, and Sutirth Dey. “A New Complexity Measure for Time Series Analysis and Classification.” The European Physical Journal Special Topics 222, no. 3–4 (July 2013): 847–60. https://doi.org/10.1140/epjst/e2013-01888-9.

## How to use
The simplest way right now is to clone this repository and use inside a conda or a pip + virtualenv environment. The only requirements for proper functionality of the entire package are `numpy` and `cython`. After cloning many functions implemented in Cython need to be compiled.

For running tests (strongly recommended), additional packages need to be installed using pip.

### Installation
0. Create a fresh conda or pip-based environment with `numpy` and `cython` packages (Skip if already available)
1. Clone the repository and enter the local directory
```
git clone https://github.com/pranaysy/ETCPy.git
cd ETCPy
```
2. Activate conda/pip environment and compile
```
python ETC/setup.py build_ext --inplace
```
3.
*Dependencies:* None! Implemented in pure Python :)
Currently unpackaged :(

Clone the repository; should be enough.

# Demo
Please check out `demo.py` to see ETC in action.

# TODO
 - Add comments and documentation
 - Extend 1-dimensional ETC to n-dimensional ETC
 - Add functions to estimate various information theoretic measures
 - Incorporate CCC
 - Visualizations
 - Functions to work with string data from text/genome sequence/etc
 - Packaging
 - Add tests

# License
Copyright 2020 Pranay S. Yadav and Nithin Nagaraj

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
