# diversity
a simple and quick approach for calculating genetic diversity
TOBE disregarded

# NOTE
# Comments regarding code functionality
The original functioning was envisioned to be with processed illumina reads that had been aligned
and sorted and analysed using other bioinformatic pipelines. Several of the functionalities regarding
this are not worked out. Instead the functioning of this is limited to the zoo_diversity_lib.py python
module, this provides a minimal working example that can function with some incredibly specific file 
formats

# Comments regarding intention behind the code
This code was never anticipated to be opened up to the wider public, development is not continuing,
it will never be fully functional.

# Comments regarding original motivation and design
The motivation for this code was witnessing multiple (and less intensive) methods
for calculating pairwise differences (divergences) that effectively reiterated a comparison between
different (python "slices") of genetic sequences. Given the nature of the measures the prior and 
proceeding states are irrelevant to the calculations, only the changes in state are required.
Consequently the approach envisioned here is a linked-list style sliding window where every 
movement of the list results in recording the changed states at either end of the list, this
is then used for the calculation of diversity and can be used to derive other parameters.
