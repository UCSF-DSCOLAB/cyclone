# <ins>Cy</ins>tometry <ins>CL</ins>ustering <ins>O</ins>ptimization a<ins>N</ins>d <ins>E</ins>valuation (Cyclone)

This repo houses the UCSF Data Science CoLab's CyTOF pipeline which performs umap calculation and clustering, and also outputs many QC and annotation focused visualizations.

Our team: Ravi Patel, Bushra Samad, Rebecca Jaszczak, Dan Bunis, Tristan Courau, Nayvin Chew, Chris Im, Arjun Rao, Lia Avesanyan, Lenny lupin, Alexis Combes and Gabi Fragiadakis

The pipeline is run primarily via the `cytof_pipeline.R` script, but we have also put together the `cyclone` R package to provide access to documentation and to help users get all of the pipeline's package dependencies installed.

---

**General Information & Instructions for Installation and Usage:** These are maintained within `cyclone`. See the [package vignette](cyclone/vignettes/Running.md) to get started.

**Citation Instructions**: We are working to put together a manuscript describing the pipeline, but it is not yet ready. For now, please reach out to us for up-to-date citation instructions if you have used the pipeline for work that you plan to publish.

**Have questions or want to contribute?**: You can reach us and/or contribute to this work by either 1) opening up an 'Issue' above to describe any bugs or feature requests in a public spot where other users might see, or 2) reaching out to a development team member by email. See the cyclone vignette for contact details. Or, 3) if you have code or documentation updates which you would like us to consider implementing directly, refer to the "Direct Contribution Instructions" below to contribute via a 'Pull Request'.

---

#### Direct Contribution Instructions
0. (Optional but recommended) Create an issue describing your proposed changes to run your idea by our team.
1. Fork this repo.
2. For contributions within the cyclone package itself, open the Cyclone/Cyclone.Rproj file.  This should open up RStudio into the Cyclone project which will have build tools set up!  You can also open RStudio and tell it to open the project by navigating to that file.  Also see the [cyclone README](cyclone/README.md) for some extra package dev instructions.
3. Commit edits to your fork of the repo.
4. Create a `Pull Request` into the `main` branch and describe the PR's goals and changes in either the PR's description or a first comment.
