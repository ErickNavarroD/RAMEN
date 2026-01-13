# Contributing to RAMEN

Hi there! This outlines how to propose a change to RAMEN. First of all,
thanks for considering contributing to our package! It’s people like you
that make it rewarding for us - the project maintainers - to work on
RAMEN. 😊

For a detailed discussion on contributing to this and other tidyverse
packages, please see the [development contributing
guide](https://rstd.io/tidy-contrib) and our [code review
principles](https://code-review.tidyverse.org/).

There are many ways you can contribute to this project (see the [Open
Source Guide](https://opensource.guide/how-to-contribute/)). Here are
some of them:

## Engage with the package

### Share the ideas

Think RAMEN is useful? Let others discover it, by telling them in
person, via BlueSky or a blog post.

Using RAMEN for a paper you are writing? Consider [citing
it](https://link.springer.com/article/10.1186/s13059-025-03864-4).

### Ask a question

Using RAMEN and got stuck? Browse the \[documentation\]\[website\] to
see if you can find a solution. Still stuck? Post your question as an
\[issue on GitHub\]\[new_issue\]. While we cannot offer user support,
we’ll try to do our best to address it, as questions often lead to
better documentation or the discovery of bugs.

Want to ask a question in private? Contact the package maintainer by
\[email\]\[[mailto:email](mailto:email)\].

### Propose an idea 💡

Have an idea for a new our_package feature? Take a look at the
\[documentation\]\[website\] and \[issue list\]\[issues\] to see if it
isn’t included or suggested yet. If not, suggest your idea as an \[issue
on GitHub\]\[new_issue\]. While we can’t promise to implement your idea,
it helps to:

- Explain in detail how it would work.
- Keep the scope as narrow as possible.

See below if you want to contribute code for your idea as well.

## Improve the documentation

Noticed a typo on the website? Think a function could use a better
example? Good documentation makes all the difference, so your help to
improve it is very welcome!

You can fix typos, spelling mistakes, or grammatical errors in the
documentation directly using the GitHub web interface, as long as the
changes are made in the *source* file. This generally means you’ll need
to edit [roxygen2
comments](https://roxygen2.r-lib.org/articles/roxygen2.html) in an `.R`,
not a `.Rd` file. You can find the `.R` file that generates the `.Rd` by
reading the comment in the first line.

## Bigger changes

If you want to make a bigger change, it’s a good idea to first file an
issue and make sure someone from the team agrees that it’s needed.

### Report a bug

Using our_package and discovered a bug? That’s annoying! Don’t let
others have the same experience and report it as well in an \[issue on
GitHub\]\[new_issue\] so we can fix it. If you’ve found a bug, please
file an issue that illustrates the bug with a minimal
[reprex](https://www.tidyverse.org/help/#reprex) (this will also help
you write a unit test, if needed). See our guide on [how to create a
great issue](https://code-review.tidyverse.org/issues/) for more advice.
Please provide as well your operating system name and version (e.g. Mac
OS 10.13.6), and any details about your local setup that might be
helpful in troubleshooting.

### Pull request process

We try to follow the [GitHub
flow](https://guides.github.com/introduction/flow/) for development.

- Fork the package and clone onto your computer. If you haven’t done
  this before, we recommend using
  `usethis::create_from_github("ErickNavarroD/RAMEN", fork = TRUE)`.

- Install all development dependencies with
  `devtools::install_dev_deps()`, and then make sure the package passes
  R CMD check by running `devtools::check()`. If R CMD check doesn’t
  pass cleanly, it’s a good idea to ask for help before continuing.

- Create a Git branch for your pull request (PR). We recommend using
  `usethis::pr_init("brief-description-of-change")`.

- Make your changes, commit to git, and then create a PR by running
  `usethis::pr_push()`, and following the prompts in your browser. The
  title of your PR should briefly describe the change. The body of your
  PR should contain `Fixes #issue-number`.

- For user-facing changes, add a bullet to the top of `NEWS.md`
  (i.e. just below the first header). Follow the style described in
  <https://style.tidyverse.org/news.html>.

### Code style

- New code should follow the tidyverse [style
  guide](https://style.tidyverse.org). You can use the
  [styler](https://CRAN.R-project.org/package=styler) package to apply
  these styles, but please don’t restyle code that has nothing to do
  with your PR.

- We use [roxygen2](https://cran.r-project.org/package=roxygen2), with
  [Markdown
  syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html),
  for documentation.

- We use [testthat](https://cran.r-project.org/package=testthat) for
  unit tests. Contributions with test cases included are easier to
  accept.

## Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.
