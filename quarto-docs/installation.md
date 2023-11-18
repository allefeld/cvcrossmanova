---
title: Installation
subtitle: How to download and install the toolbox and its dependencies
order: 1
---

The CvCrossManova toolbox has been developed with Matlab R2023b, but should work with later and slightly earlier versions; the most recently introduced feature used by the code is from R2021a.

The core functionality does not depend on other toolboxes. However, functions which read data modeled in SPM12 (`ccmRegion`, `ccmSearchlight`, `ModeledData.fromSPM`) depend on it being installed. You can download SPM12 from the pages of the [Functional Imaging Lab at UCL](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/). The toolbox was developed using r7771, but other versions of SPM12 should work, too. It has not been tested with SPM8 or earlier.

To install the toolbox, download the file `cvcrossmanova-v#.#.#.zip` from the latest GitHub [release](https://github.com/allefeld/cvcrossmanova/releases). Unzip the file into a directory, and make sure that the created directory `cvcrossmanova-v#.#.#/` is on the [Matlab search path](https://www.mathworks.com/help/matlab/matlab_env/what-is-the-matlab-search-path.html). If you install a newer version, make sure that the old version is no longer on the path.

The code for the toolbox as well as this documentation is maintained in a git [repository](https://github.com/allefeld/cvcrossmanova) on GitHub. If you need the latest not yet released version, [clone](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository) the repository with
```bash
git clone https://github.com/allefeld/cvcrossmanova.git
```
and copy the files in the `toolbox/` subdirectory into a directory on the Matlab search path.

If you found a bug or have an idea for improvement, you can create an [issue](https://github.com/allefeld/cvcrossmanova/issues) on the repository.