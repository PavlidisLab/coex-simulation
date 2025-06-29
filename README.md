# A theoretical study on the detection of dynamic gene regulation in bulk tissue transcriptomes.
C. Pan Chu, Alexander Morin, Paul Pavlidis

Thousands of studies have used co-expression analysis of bulk tissue samples to probe gene regulation. However, the extent that intracellular regulatory signals are present in these data is unclear. Specifically, we lack clarity of the factors that promote or impede the propagation of regulatory signals from the single cell level to the bulk tissue level. To bring these issues into focus, we developed a novel computational simulator, grounded in real data, to explore the theoretical relationship between events in single cells and bulk tissue expression profiles, and clarify the conditions required for the propagation of intracellular regulatory signals in complex tissues such as the brain. Our simulator first generates single cell expression profiles and subsequently samples and aggregates these single cells to produce bulk tissue expression profiles. Using this framework, we found that there are very specific and unlikely conditions under which intracellular dynamic regulatory signals can be propagated to the bulk tissue level. For the most part, such regulatory relationships, however strong at the single cell level, are unlikely to be detectable. Our results provide a quantitative explanation for why regulatory network inference from co-expression has proved challenging - even with the assistance of other data modalities - and gives the scientific community a set of tools to further explore these issues in both single-cell and bulk tissue data.

# Gene Expression Simulator

The simulator code is contained in the_simulator.R. 

This R-based simulator generates gene expression data. It incorporates **subject-level variability**, **cell-level heterogeneity**, and **gene co-expression patterns** to produce synthetic datasets.

## Overview

The simulator learns statistical models from user-provided reference data to generate new, synthetic expression matrices. Its core capabilities, accessed via an `api` object, allow users to:

*   **Initialize** the simulator.
*   **Provide reference datasets** for marginal distributions and cell type expression profiles.
*   **Define gene co-expression programs** at subject and cell levels.
*   **Fit statistical models**, specifically Gamma distributions, to capture gene-wise mean and variance relationships from reference data.
*   **Simulate new gene expression data** by combining baseline cell simulations with subject-specific expression biases and co-expression patterns.

## How to Use the Simulator to Generate New Data

Generating new gene expression data with this simulator involves a sequence of steps using functions exposed through the `api` list.

### 1. Initialize the Simulator Object

Start by creating a simulator instance.

```R
# Call initSimulator to get the 'this' object and the 'api' functions
simulator_environment <- initSimulator()
simulator_obj <- simulator_environment$this
api <- simulator_environment$api

2. Set Reference Data
The simulator learns parameters from provided reference datasets.
    •Reference Expression Matrices (exprmat): Used to estimate marginal distributions.
        ◦Subject-level (sbj): For subject-level expression counts.
        ◦Cell-level (cel): For individual cell expression counts.
    •Cell Type Expression Profiles (cteprf): Defines mean expression per gene across cell types.
    •Co-expression Programs (coexPrograms): Specifies gene groups that co-express and their correlation.
        ◦Subject-level (sbj): For co-expression across subjects. Cell types sharing a program synchronize at the subject level.
        ◦Cell-level (cel): For co-expression within individual cells.

3. Fit Marginal Distribution Models
This step fits Gamma distribution models for each gene. A third-degree polynomial model (lm(variance ~ poly(mean, 3, raw = TRUE))) predicts variance from mean; residuals add variability.
simulator_obj <- api$fitMds(simulator_obj)

4. Simulate Baseline Cells (Cell-Level Variability)
This generates initial cell-level expression, incorporating cell-level co-expression. Internal private$utils$generateVals samples from Gamma distributions with co-expression.
# nSubject: number of subjects; nCell: number of cells per subject per cell type
simCells_output <- api$simBseLnCels(simulator_obj, nSubject = 10, nCell = 100)

5. Simulate Subject-Level Means (Subject-Level Variability)
This generates subject-level mean expression for each gene and cell type, including subject-level co-expression patterns. Internal private$utils$generateVals is used for sampling.
# nSubject: number of subjects
simSbjs_output <- api$simSbjLvMeans(simulator_obj, nSubject = 10)

6. Compute Combined Cell-Level Parameters
Subject-level means adjust cell-level distribution models, creating unique cell-level Gamma distribution parameters for each gene, cell type, and subject. The original cell-level distribution's mean is adjusted while maintaining the relative variance (coefficient of variation, variance / mean) for each gene.
cel_params_combined <- api$computeCelParams(simulator_obj, sbjLvMeans = simSbjs_output$simSbjs)

7. Convert Cell-Level Distributions to Final Expressions
Baseline cells are transformed to reflect subject-level variability. This involves converting baseline cell expressions to p-values using their original Gamma distribution parameters, then transforming these p-values back into new expression values using the subject-specific cell-level parameters.
final_exprmat_output <- api$convertCelLvDist(
    simCells = simCells_output$simCells,
    celMdParamsOrig = simulator_obj$mdParams$cel,
    celMdParamsNew = cel_params_combined
)

Additional Utility Function
    •api$GENERATE_CC_SPECS(nBkSamples, nTotalCells, baselineProps, sdFrac): Generates cellular composition specifications. It produces a matrix specifying cell counts per cell type for subjects, using baselineProps (proportion of each cell type) and sdFrac (standard deviation fraction) for sampling from a normal distribution.
# Example: Generate specs for 5 subjects, 1000 total cells, with baseline proportions
# baseline_props <- c(typeA = 0.5, typeB = 0.5)
# cell_comp_specs <- api$GENERATE_CC_SPECS(nBkSamples = 5, nTotalCells = 1000, baselineProps = baseline_props, sdFrac = 0.1)
