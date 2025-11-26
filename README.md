# CathodeDataExtractor

------------

__astrochemreactionextractor__  is a lightweight document-level information extraction pipeline that can automatically extract chemical reaction equations from chemical literature, along with the names of species related to the reactions and kinetic parameters.

Prerequisites
------------

- __git__: Install `git` for downloading the package https://git-scm.com/book/en/v2/Getting-Started-Installing-Git.
- __conda__: Dependency manager used to create Python virtual environments https://conda.io/projects/conda/en/latest/user-guide/install/index.html.

## Installation

------------
Download and go to the directory of the repository.
```
git clone https://github.com/Dingyun-Huang/chemdataextractorTADF.git
cd chemdataextractorTADF
```

Create and activate a new Python 3.11 environment.
```
conda create --name my_env python=3.11
conda activate my_env
```

When you are in the repository directory, install astrochemreactionextractor.
```
pip install astrochemreactionextractor
```

## Quick start

------------
#### Extract from documents

```python
from astrochemreactionextractor.reaction_extraction_pipe import ReactionPipeline

input_path = './papers_for_extraction'
output_path = './output'

pipline  = ReactionPipeline(input_path, output_path)
pipline.extract()
```
> 

## Issues?

------------
You can either report an issue on GitHub or contact me directly. Try buhl@zhejianglab.org.

## Citing 

------------
If the source code turns out to be helpful to your research, please cite the following work:
