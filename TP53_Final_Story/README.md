# TP53 Final Story
Using the CPTAC package, we analyzed data from endometrial cancer patients in order to determine the effects of TP53 mutation on protein expression and phosphorylation in endometrial cancer cells. We focused our analyses on comparing three groups: cancer wildtype (no TP53 mutation), TP53 mutation within a mutational "hotspot" (the DNA binding domain of TP53), and other TP53 mutations. The contents of the individual Jupyter notebooks are as follows:

(1) <b>TP53_Cis:</b> Our first step was to determine the cis effects of TP53 mutation on p53.

(2) <b>TP53_Interacting_List:</b> We compiled a list of proteins known to interact with TP53 (according to Uniprot and String) and tested each one for a significant change in protein abundance in the presence of TP53 mutation. We also analyzed if the location of TP53 mutation (hotspot/non-hotspot) made a difference.

(3) <b>TP53_All_Proteins:</b> We decided to broaden our analysis and look at the potential effects of TP53 mutation on all of the proteins in our dataset (~11,000) using methods very similar to those used in (2). This revealed that TP53 mutation had a significant impact on many proteins that were not previously associated with p53.

(4) <b>TP53_Interesting_Protein_Analysis:</b> A more in-depth proteomic and phosphoproteomic analysis of interesting proteins found in notebooks (2) and (3), namely AURKA, BORA, XPO1, CDK1, TPX2, PLK1, CABLES1, and XRN2.

(5) <b>TP53_Clinical:</b> An analysis of the association between TP53 mutation and various clinical outcomes.

(6) <b>TP53_Transcriptomics:</b> A brief look into the analyses we did of the impact on TP53 mutation on RNA expression of various proteins of interest. This is more of a side note, supplementary to our main analyses but not covered in the manuscript.
