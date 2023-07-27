# <ins>Cy</ins>tometry <ins>CL</ins>ustering <ins>O</ins>ptimization a<ins>N</ins>d <ins>E</ins>valuation (Cyclone)

This repo houses the UCSF Data Science CoLab's CyTOF pipeline which performs umap calculation and clustering, and also outputs many QC and annotation focused visualizations.

Our team: Ravi Patel,  Rebecca Jaszczak, Chris Im, Nicholas Carey, Tristan Courau, Dan Bunis, Bushra Samad, Lia Avanesyan, Nayvin Chew, Sarah Stenske, Jillian M. Jespersen, Jean Publicover, Austin Edwards, Mohammad Naser, Arjun Rao, Lenny Lupin-Jimenez, Matthew Krummel, Stewart Cooper, Jody Baron, Alexis Combes, and Gabriela Fragiadakis

The pipeline is run primarily via the `cytof_pipeline.R` script, but we have also put together the `cyclone` R package to provide access to documentation and to help users get all of the pipeline's package dependencies installed.

---

**General Information & Instructions for Installation and Usage:** These are maintained within `cyclone`, but we've added knit'd versions of all vignettes to this repo. See the ["Running" vignette](vignettes/Running.md) to get started, and the ["FollowUp" vignette](vignettes/FollowUp.html) for ideas of what can be done with cyclone outputs after running the pipeline.

**Citation Instructions**: We have put together a manuscript describing the pipeline.  [bioRxiv preprint is here]( https://www.biorxiv.org/content/10.1101/2023.03.08.531782v1), peer-reviewed version is coming soon!

**Have questions or want to contribute?**: You can reach us and/or contribute to this work by either 1) opening up an 'Issue' above to describe any bugs or feature requests in a public spot where other users might see, or 2) reaching out to a development team member by email. See the cyclone vignette for contact details. Or, 3) if you have code or documentation updates which you would like us to consider implementing directly, refer to the "Direct Contribution Instructions" below to contribute via a 'Pull Request'.

---

#### Direct Contribution Instructions
0. (Optional but recommended) Create an issue describing your proposed changes to run your idea by our team.
1. Fork this repo.
2. For contributions within the cyclone package itself, open the Cyclone/Cyclone.Rproj file.  This should open up RStudio into the Cyclone project which will have build tools set up!  You can also open RStudio and tell it to open the project by navigating to that file.  Also see the [cyclone README](cyclone/README.md) for some extra package dev instructions.
3. Commit edits to your fork of the repo.
4. Create a `Pull Request` into the `main` branch and describe the PR's goals and changes in either the PR's description or a first comment.
