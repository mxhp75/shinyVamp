## Removing miRNA from haemolysis metric calculation

***

Prior to calculation of the haemolysis metric, the difference between the geometric mean of signature and background miRNA,
we must first remove from the siganture group any miRNA that are _a priori_ understood to be different between any group of interest.

***

Briefly, if in the context of the experiment curently being investigated a differential expression has been performed, any miRNA that are statistically differentially expressed between the groups and that interestect with the 20 miRNA in the signature set should be selected using the radio buttons. Selecting these miRNA will remove them from  the signature set and subsequent calculations. If in there has been no differential expression analysis but there are miRNA known from the literature to be differential between the groups, these miRNA should also be selected, removing them from future calculations.

***

For example, if a differential expression has been performed between hypothetical _groupA_ (control) and _groupB_ (treatment) and **hsa-miR-106b-3p**, **hsa-miR-9-5p** and **hsa-miR-186-5p** are statistically differentially expressed between groupA and groupB, select the radio button for **hsa-miR-106b-3p** and **hsa-miR-186-5p**. This step helps to ensure any issues with haemolysis detection are not confounded with the groups of interest.

