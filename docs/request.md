# Revised analysis request

Compare haplotypes across Diabetes, Prediabetes and Control using the existing project data and paper-inspired Figure C as the visual reference.

1. Reproduce Figure C with accurate labels for the quantities it actually measures.
2. Add a comparison figure showing Diabetes versus Control and Prediabetes versus Control. Show eligible haplotype carrier differences, confidence intervals and adjusted p-values; label unavailable comparisons explicitly.
3. Draw a separate regional marker picture for every GWAS-selected region, using all available qualifying regional markers rather than only the lead marker or a single example region. Provide an index and a separate interpretation for each region.
4. Report eligible overall and pairwise haplotype comparisons, effect sizes, confidence intervals and multiple-testing correction, distinguishing descriptive chromosome-copy frequencies from participant-level inference.
5. Investigate candidate causal variants using available genetic and external evidence. Distinguish association from causality and explain what evidence is missing when localization is unresolved.
6. Provide complete reproducible Python and independent R implementations, documentation and AGENTS.md in both GitHub repositories.
7. Publish code and generic methods only in public `Homabazireh/fine-mapping`. Store all project data, results, figures and cohort interpretations only in private `Homabazireh/fine-mapping-2`, with verified archive checksums and restoration instructions. Exclude credentials, environments and caches. Do not invent inputs, results or causal conclusions.
