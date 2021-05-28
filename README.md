# diversity
a simple and quick approach for calculating genetic diversity

# NOTE
this is currently none functional and takes code from an intermediate stage of development
initially this code was written for more or less complete sequences. This version is 
being developed to be used with short reads that have been aligned using a reference.
The basic input is envisioned to be FASTA sequences with the names conveying the start 
and end of the short sequences in relation to the position on the reference. One of the basic
concepts in this is that errors are expected to be only recorded singly and can therefore
be isolated simply.

## NOTE
this code is in a rough state and was never anticipated to be opened up to the wider public,
it is also currently not operational. Earlier versions were used for preliminary analyses,
however this was not with short reads but using multiple references.

## NOTE
the motivation for this code was witnessing multiple (and less programmer intensive) methods
for calculating pairwise differences (divergences) that effectively reiterated a comparison between
different (python "slices") of genetic sequences. Given the nature of the measures the prior and 
proceeding states are irrelevant to the calculations, only the changes in state are required.
Consequently the approach envisioned here is a linked-list style sliding window where every 
movement of the list results in recording the changed states at either end of the list, this
is then used for the current calculation of diversity.
