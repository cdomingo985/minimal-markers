MINIMUM MARKER SETS 

min_marker_set.R

Program purpose: 
    
    From a SNP large dataset genotypied collection of varieties a list of minimum sets 
    of markers that could unequivocally identify each of the varities 

SNP markers have provided a reliable and stable tool for plant genotyping for plant variety identification and protection. The availability of small and low-cost SNP panels to accelerate the identification of the cultivated rice varieties should be beneficial for breeders, seed certification entities and rice industry. With the intention of providing of such a facility, we first developed a simple and easy-handle bioinformatics tool based on the widely used and freely available software R to generate small sets of SNPs that can discriminate varieties, by selecting markers from a larger genotyping dataset. The provided R-based algorithm can be applied to other genomic resources to develop core sets of discriminating markers. 

Using genotypic dataset to discriminate the varieties, first SNPs that showed deletions in some accession are removed. Second, it identifies, if any, pairs or groups of varieties that cannot be discriminated by these markers. If such pair of varieties was found, they were treated as a single accession. At this step, the algorithm generated a list of these pairs of non-differentiable accessions. In the third step it is performed pairwise comparisons among the genotypes to find the group of markers that discriminate each pair. For these markers, the algorithm calculates those that differentiate the same pairs of varieties, that is, the redundant markers that are removed. One marker, the first in the list, is selected to proceed to the following step and the remaining are removed for calculations. Finally, pairwise comparisons were performed by iteration, which produced different combinations of SNP markers that distinguished the varieties in a collection 


The use of this algorithm has been published in PONE-S-23-05902

Authors:    
    Reig Valiente, Juan Luis
    Domingo Velasco, Concha
    Garcia-Romeral, Julia
