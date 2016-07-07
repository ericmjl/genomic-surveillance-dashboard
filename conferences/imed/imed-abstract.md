---
title: A real-time surveillance dashboard for monitoring viral phenotype from sequence.
author:
- name: Eric J. Ma
  affiliation: MIT
- name: Jonathan A. Runstadler
  affiliation: MIT
---

# Purpose

Genome sequencing has become routine in outbreak surveillance. Apart from its utility in tracking evolutionary trajectories, genome sequences also contain information necessary for predicting biochemical phenotypes of outbreak pathogens, such as drug resistance and antigenic distance from vaccine strains. Having the capability to predict and visualize a pathogen's phenotype from its genome, in real-tiem, and at a location near the source of the outbreak, would provide epidemiologists with additional useful data that could help tailor outbreak responses and appropriately mobilize limited resources.

# Materials and Methods

We provide a design for an extensible digital dashboard for infectious disease genomic surveillance. This open-source dashboard is written in the Python programming language. The application back-end is written using the Flask web framework. Automated machine learning model and parameter selection is enabled using the `numpy` and `scikit-learn` packages. Data visualization is provided by the `bokeh` package. Data are stored in SQLite databases. Automated tests using `py.test` are used to test for data and application integrity. The dashboard's backend is modularly designed, with separate "microservices" for the web interface, data ingestion and preprocessing, machine learning, and visualization. This design enables the addition of new pathogens, phenotypes, machine learning models, and visualization types as needed.

# Results

We have developed a proof-of-concept dashboard, using sequence-phenotype data from the HIV drug resistance database and the Los Alamos HIV Sequence Database.The dashboard takes in a new virus' sequence, and returns a visualization of its predicted drug resistance profile. We show how research groups can programmatically submit standardized phenotype data to update the database and machine learning models. We demonstrate the portability of the data and models. We also provide examples on how to extend the dashboard.

# Conclusion

New "gold-standard" genotype-phenotype data and integrated prediction systems are required to fully realize the utility of genomic sequencing data by connecting genotype to phenotype. We have designed this dashboard, and implemented a proof-of-concept, as part of our efforts to realize this vision.
