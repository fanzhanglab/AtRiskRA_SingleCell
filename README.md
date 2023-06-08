# Deep immunophenotyping reveals circulating activated lymphocytes in individuals at risk for rheumatoid arthritis

Rheumatoid arthritis (RA) is a systemic autoimmune disease with currently no effective prevention strategies. Identifying pathogenic immune phenotypes in ‘At-Risk’ populations prior to clinical disease is crucial to establishing prevention strategies. Here, we applied mass cytometry to deeply characterize the immunophenotypes in blood from At-Risk individuals of clinical subpopulations based on family history and autoantibody status (n=52), established RA (n=67), and healthy controls (n=48). Through integrative and disease association analyses, we uncovered significant cell expansions in At-Risk individuals compared with controls, including CCR2+ T helper cells, T peripheral helper cells, Type 1 T helper cells, CXCR5+CD8+ T cells. We further validated the At-Risk associations of T cell phenotypes using our validation cohort with 57 At-Risk and 23 healthy individuals. In addition, we found that CD15+ classical monocytes were especially expanded in ACPA-negative At-Risk, and an activated PAX5low naïve B cell population expanded in ACPA-positive individuals who also had an FDR with RA. Moreover, we developed a “RA immunophenotype score” classification method based on the degree of enrichment and abundance of cell states relevant to established RA using mixed effect modeling and logistic regression, and demonstrated that this score significantly distinguished At-Risk individuals from the control (p=0.039 and AUC>0.6). In all, we systematically characterized immunophenotypical differences among At-Risk subpopulations and their abundance comparing with established RA using single-cell proteomics; Further, our classification model may provide a promising approach towards uncovering the unknown pathogenesis of RA with the goal to develop preventive strategies.

**Reference: Inamo J et al. Deep immunophenotyping reveals circulating activated lymphocytes in individuals at risk for rheumatoid arthritis. ([preprint](https://XXX))**


# Study design
![image](./images/CyTOF_workflow.png)

# Highlights
- We identified immune cell phenotypes relevant to the preclinical phase of RA
- We developed "RA immunophenotype score" schema to estimate the degree of RA-relevant immunophenotypes for individuals
- We provide large-scale reference dataset for single-cell proteomics analysis

# Results summary
![image](./images/results.jpg)
Identifications of specific T cell populations that were associated with At-Risk. A. Cells in UMAP are colored in red (expansion) or blue (depletion) and p-value is shown as well. B. Distributions of cell neighborhood correlations and odds ratio. Error bars represent 95% confidence intervals. 

![image](./images/RA_immunophenotype_score.jpg)
Classifications for At-Risk individuals and established RA individuals. A. RA immunophenotype score utilizing RA-specific cell type abundances to quantify and distinguish At-Risk individuals from control. For each cell type, all p-values from the covarying neighborhood analysis test were p = 1e-3. We incorporated clusters that are significantly associated with RA (adjusted p < 0.05) to model the RA immunophenotype score. We calculated RA immunophenotype score based on cell type abundances multiplied by corresponding major cell type proportions and enrichment scores for each cell type, B. Distribution of RA immunophenotype score across individual samples from RA, At-Risk, and controls; **** p < 0.0001, * p < 0.05, C. Receiver operating characteristic (ROC) analysis to evaluate the classification performance of RA immunophenotype score in distinguishing At-Risk from control. Areas under the curve (AUC) with 95% confidence intervals were described. All the analyses are adjusted for age and sex.