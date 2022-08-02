# CyTOF CLustering Optimization aNd Evaluation (Cyclone)

This repo houses the UCSF Data Science CoLab's CyTOF pipeline which performs umap calculation and clustering, and also outputs many QC and annotation focused visualizations.

Our team: Ravi Patel, Bushra Samad, Rebecca Jaszczak, Dan Bunis, Tristan Courau, Nayvin Chew, Chris Im, Arjun Rao, Lia Avesanyan, Lenny lupin, Alexis Combes and Gabi Fragiadakis

The pipeline is run primarily via the `cytof_pipeline.R` script, but we have also put together the `cyclone` R package to provide access to documentation and to help users get all of the pipeline's package dependencies installed.

---

**General Information & Instructions for Installation and Usage:** These are maintained within `cyclone`. See the [package vignette](vignettes/Running.md) to get started.

**Citation Instructions**: We are working to put together a manuscript describing the pipeline, but it is not yet ready. For now, please reach out to us for up-to-date citation instructions if you have used the pipeline for work that you plan to publish.

**Have questions or want to contribute?**: You can reach us and/or contribute to this work by either 1) opening up an 'Issue' above to describe any bugs or feature requests in a public spot where other users might see, or 2) reaching out to a development team member by email. See the cyclone vignette for contact details. Or, 3) if you have code or documentation updates which you would like us to consider implementing directly, refer to the "Direct Contribution Instructions" below to contribute via a 'Pull Request'.

---

# To help develop this package with RStudio:

Instructions here, by Dan, are for anyone who wants to help develop the `cyclone` package.

0. (Optional but recommended) Create an issue describing your proposed changes to run your idea by our team.
1. Fork this repo.
2. For contributions within the cyclone package itself, open the Cyclone/Cyclone.Rproj file.  This should open up RStudio into the Cyclone project which will have build tools set up!  You can also open RStudio and tell it to open the project by navigating to that file.
3. Commit edits to your fork of the repo.
4. Create a `Pull Request` into the `main` branch and describe the PR's goals and changes in either the PR's description or a first comment.

## Writing Documentation:
I've always used roxygen for documentation management.
It allows for code-based documentation creation which I've found pretty easy to use.
Roxygen documentation gets written right into a package's `.R` files in the `R/` directory, with any lines starting with `#'` being the roxygen lines.
For a normal package, that let's you document functions in the same spot that you define them.
Of course, we won't always be doing that exactly here... You can also follow documentation chunks with a `NULL` in order to reference code that is maintained elsewhere.
**I got a few started already [here](R/checkpoints.R)!**

You can check out [the dittoSeq repo](https://github.com/dtm2451/dittoSeq/tree/master/R) for some of my other examples, or the roxygen2 vignettes [here](https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html).

### Building the documentation:
How the project is set up, documentation is re-built any time that the `'Build ...'` or `'Check'` actions are made (More on those below). That's what I end up relying on mostly in dittoSeq dev, but you can also use `roxygen2::roxygenise()`, `devtools::build()`, or `Ctrl + Shift + D`!

## Writing the Vignette:
It's an `.Rmd` file in the `vignettes` folder.  Rmd files give a rendered text md format surrounding code chunks that can be controlled with arguments provided to the \`\`\`{r} that they start with.

### Building the vignette:
With the `.Rmd` open, click the `Knit` button in the menu bar of the editor panel.  This will create a local copy of the vignette which you can view.  (That file does not need to be pushed to the repo. The vignette will be automatically constructed within any `'Build ...'` action and included in the packaged output.)   

## Building / Checking / Installation
Build tools are in the top right panel of RStudio, in the `Build` tab.

The 3 actions I use most often:
- Build a source version of the package with `More > Build Source Package`
- Install from the source version which ^^ will build with `Install and Restart`
- Run Checks: The `Check` button performs lots of tests and is probably what I've relied on most of anything.

To view the documentation yourself, Build, then install, then `?<target_docs>`

To run automated checks which include checks of the documentation, use the Check feature.

**I learned by doing and since this is a git repo, we can always undo if something breaks!**
