# diversity
a simple and quick approach for calculating genetic diversity

# NOTE
# Comments regarding code functionality
this is currently non-functional and takes code from an intermediate stage of development,
initially this code was written for more or less complete sequences. In contrast, this version 
was developed to be used with short reads aligned to a reference. The basic input was 
envisioned to be FASTA sequences with the names that conveyed the start and end of the 
sequences in relation to the position on the reference. One of the basic concepts of this 
approach was that errors are expected to occur singly within the reads and can therefore
be easily isolated.

# Comments regarding intention behind the code
this code was never anticipated to be opened up to the wider public, development is not continuing
and it is not envisioned to be functional. The first versions (note the previous comment) were used 
for preliminary analyses.

# Comments regarding original motivation and design
the motivation for this code was witnessing multiple (and less intensive) methods
for calculating pairwise differences (divergences) that effectively reiterated a comparison between
different (python "slices") of genetic sequences. Given the nature of the measures the prior and 
proceeding states are irrelevant to the calculations, only the changes in state are required.
Consequently the approach envisioned here is a linked-list style sliding window where every 
movement of the list results in recording the changed states at either end of the list, this
is then used for the calculation of diversity.
