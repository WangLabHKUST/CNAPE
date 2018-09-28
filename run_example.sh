#! /bin/bash

#the cnape.R takes two arguments:
#the first one is your expression matrix
#the format is examplified in the example/example.geneExpression.RPKM.txt
#the RPKM values should be prepared following the TCGA pipeline

#the second argument is the prefix of the output.
#Two files will be written to the output:
#${prefix}.chromosome_level.cna.txt       ${prefix}.arm_level.cna.txt
Rscript cnape.R example/example.geneExpression.RPKM.txt example/example
