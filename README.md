# ETCPy
Compute Effort To Compress (ETC) using the NSRPS algorithm in Python.  
  
Currently computes the exact ETC on a 1-dimensional symbolic sequence. Convenience functions provided for converting data to symbolic sequence.  
  
***Reference paper:***  
Nagaraj, Nithin, Karthi Balasubramanian, and Sutirth Dey. *A New Complexity Measure for Time Series Analysis and Classification.* The European Physical Journal Special Topics 222, no. 3–4 (July 2013): 847–60. https://doi.org/10.1140/epjst/e2013-01888-9.
  
# Installation
Dependencies: `Numba` and `NumPy`
Currently unpackaged. Clone the repository, install NumPy and Numba using conda.
  
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
