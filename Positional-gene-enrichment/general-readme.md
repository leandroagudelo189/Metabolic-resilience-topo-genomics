### The code in this repository has been updated and modified based on the methodology from the following study, as the hosting of their web user interface has been discontinued. Full details of the methodology can be found in the referenced publication.


<br>
<br>
<br>

[Nucleic Acids Res.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2367735/#) 2008 Apr; 36(7): e43.  
Published online 2008 Mar 16\. doi: [10.1093/nar/gkn114](https://doi.org/10.1093%2Fnar%2Fgkn114)

# Positional gene enrichment analysis of gene sets for high-resolution identification of overrepresented chromosomal regions

[Katleen De Preter](https://pubmed.ncbi.nlm.nih.gov/?term=De%20Preter%20K%5BAuthor%5D),1,\* [Roland Barriot](https://pubmed.ncbi.nlm.nih.gov/?term=Barriot%20R%5BAuthor%5D),2 [Frank Speleman](https://pubmed.ncbi.nlm.nih.gov/?term=Speleman%20F%5BAuthor%5D),1 [Jo Vandesompele](https://pubmed.ncbi.nlm.nih.gov/?term=Vandesompele%20J%5BAuthor%5D),1 and [Yves Moreau](https://pubmed.ncbi.nlm.nih.gov/?term=Moreau%20Y%5BAuthor%5D)2

[Author information](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2367735/#) [Article notes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2367735/#) [Copyright and License information](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2367735/#) [PMC Disclaimer](https://www.ncbi.nlm.nih.gov/pmc/about/disclaimer/)

[Go to:](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2367735/#)  Abstract: The search for feature enrichment is a widely used method to characterize a set of genes. While several tools have been designed for nominal features such as Gene Ontology annotations or KEGG Pathways, very little has been proposed to tackle numerical features such as the chromosomal positions of genes. For instance, microarray studies typically generate gene lists that are differentially expressed in the sample subgroups under investigation, and when studying diseases caused by genome alterations, it is of great interest to delineate the chromosomal regions that are significantly enriched in these lists. In this article, we present a positional gene enrichment analysis method (PGE) for the identification of chromosomal regions that are significantly enriched in a given set of genes. The strength of our method relies on an original query optimization approach that allows to virtually consider all the possible chromosomal regions for enrichment, and on the multiple testing correction which discriminates truly enriched regions versus those that can occur by chance. We have developed a Web tool implementing this method applied to the human genome ([http://www.esat.kuleuven.be/\~bioiuser/pge](http://www.esat.kuleuven.be/~bioiuser/pge)). We validated PGE on published lists of differentially expressed genes. These analyses showed significant overrepresentation of known aberrant chromosomal regions.


<br>
<br>
<br>


# **Method summary for this study: Identification and functional annotation of genomic hubs by positional gene enrichments (PGE)**

To identify genomic hubs from the positional location of genes of interest, we used previously established positional gene enrichment methods (De Preter et al. 2008, NAR). Briefly, PGEs exploit topological features of gene locations to determine chromatin regions enrichments. This approach is based on calculating the hypergeometric distribution along genomic distances for gene sets. For a given region, it determines the probability of having observed genes in that region. To test statistical significance, cumulative p-value distributions are used on random simulations or false discovery rate on large gene-sets. The probability of achieving p-value enrichments that are better than chance estimates are used to define hub enrichments. The additional constraints to categorize a genomic region displaying positional enrichments for a given set of genes are as follows: having at least 2 genes; that no smaller region was found by random permutations and cumulative p value distributions; that no bigger regions with more target genes were found. This algorithm defines genomic regions by distance of query genes, followed by estimation of p-value distribution to compare chance expectation. Target genomic ranges for positional enrichment were set to 15-million base pairs (following structural relationships on interacting homotypic TADs), where PGEs with larger domains were filtered out. To obtain qualitative PGEs we used regions that were enriched within cytogenetic bands. To obtain a global algorithmic performance determined by the number of genes being queried, we additionally performed permutations with random gene sets (same number of inputs), followed by statistical comparisons using Chi-square test and p-value discrimination. For every topological gene cluster in the form of PGEs, we derived functional enrichments using Enrichr in R. FDR \<0.01 and p-value \< 0.05 was used as a threshold to select the significant enrichments for biological terms and phenotypes associated with a given region. PGEs with common biological enrichment across domain-catalogs were categorized as class fPGEs (functional PGEs hubs) for downstream analysis. These hierarchical chromatin domains were used as seeds for prioritization of structural chromatin data when indicated, and as inputs for long- and short-range structural relationships. 
