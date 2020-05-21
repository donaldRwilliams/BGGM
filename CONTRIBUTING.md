 # Contributing to BGGM development

Some of this is copied from the `ggplot2` contributing [file](https://github.com/tidyverse/ggplot2/blob/master/CONTRIBUTING.md)

# Issues
When filing an issue, the most important thing is to include a minimal reproducible example so that we can quickly verify the problem, 
and then figure out how to fix it. 

The easiest way to include data is to use dput() to generate the R code to recreate it. For example, to recreate the mtcars dataset in R, I'd perform the following steps:

Run dput(mtcars) in R
Copy the output
In my reproducible script, type mtcars <- then paste.
But even better is if you can create a data.frame() with just a handful of rows and columns that still illustrates the problem.

Spend a little bit of time ensuring that your code is easy for others to read:

make sure you've used spaces and your variable names are concise, but informative

use comments to indicate where your problem lies

do your best to remove everything that is not related to the problem.
The shorter your code is, the easier it is to understand.

You can check you have actually made a reproducible example by starting up a fresh R session and pasting your script in.

# Pull requests
To contribute a change to `BGGM`, you follow these steps:

Create a branch in git and make your changes.
Push branch to github and issue pull request (PR).
Discuss the pull request.

## Coding Conventions
Unless absolutely necessary (and when is it, really), all code should be written in base R or `c++`
