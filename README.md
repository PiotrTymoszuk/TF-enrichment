# TF-enrichment
Tools for checking transcription factor binding site enrichment in the gene set over the whole genome

The whole genome TF binding site counts were obtained with D-Light tool (https://pbwww.services.came.sbg.ac.at/?page_id=40). At the moment works only for murine genes.
The enrichment is assessed by the bootstrap technique - it's checked if there's more binding sites of the particular transcription factor in the random subsets of the whole-genome data base than in the gene set of interest.
