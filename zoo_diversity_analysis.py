#/usr/bin/env python3.7

import zoo_diversity_lib


#diversity_analyses.main(arguments.shorah_input,arguments.output_file,arguments.reference_file,arguments.input_is_directory=False,arguments.contiguous_filter=10,arguments.identity=0.85,arguments.MEGA=None)
def main(shorah_input,output_file,contiguous_filter=10,identity=0.95,MEGA=None,window_length=66,wattersons=True,coverage=0,posterior=0.95,min_haplotypes=2,shorah=False):
    fileInput = open(shorah_input,"r")
    inputSequences = zoo_diversity_lib.loadSeq(fileInput)



    alignmentA = zoo_diversity_lib.Alignment(shorah_input,window_length)
    alignmentA.set_shorah_input(shorah)
    alignmentA.set_Wattersons(wattersons)
    alignmentA.set_minimum_coverage(coverage)
    alignmentA.set_required_no_haplotypes(min_haplotypes)
    alignmentA.set_required_posterior(posterior)

    for key,val in inputSequences.items():
        alignmentA.add(key,val)

    alignmentA.start()
# the Alignment.compare(self,other,filestream) method is to be used to 
# direct the output of comparisons between the alignment objects
    alignmentA.output(output_file)

    return
