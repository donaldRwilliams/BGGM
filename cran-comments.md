# Version 2.1.5
## Addressing CRAN errors
- Adjusted autoconf.ac, dropped refence to openMP as parallelization is not implemented in the c++ file

# Version 2.1.4
## Addressing CRAN submission errors
- CRAN server running example in bggm_missing returned `Error: inv(): matrix is singular`
  Error can not be locally reproduced using R CMD check --as-cran or rhub. 
  Temporary workaround is to remove offending example

## R CMD check results

❯ checking installed package size ... NOTE
    installed size is 12.2Mb
    sub-directories of 1Mb or more:
      doc    3.3Mb
      help   1.1Mb
      libs   7.0Mb

❯ checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-Werror=format-security’ ‘-Wformat’ ‘-Wp,-D_FORTIFY_SOURCE=3’
    ‘-Wp,-D_GLIBCXX_ASSERTIONS’ ‘-march=x86-64’
    ‘-mno-omit-leaf-frame-pointer’

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

## revdepcheck results
Checked 1 rev dep: zero problems.

# Version 2.1.3
## R CMD check results

❯ checking installed package size ... NOTE
    installed size is 12.1Mb
    sub-directories of 1Mb or more:
      doc    3.4Mb
      help   1.1Mb
      libs   6.9Mb

❯ checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-Werror=format-security’ ‘-Wformat’ ‘-Wp,-D_FORTIFY_SOURCE=3’
    ‘-Wp,-D_GLIBCXX_ASSERTIONS’ ‘-march=x86-64’
    ‘-mno-omit-leaf-frame-pointer’

0 errors ✔ | 0 warnings ✔ | 2 notes ✖


# Version 2.1.2
## R CMD check results

❯ checking installed package size ... NOTE
    installed size is  8.8Mb
    sub-directories of 1Mb or more:
      help   1.1Mb
      libs   6.9Mb

❯ checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-Werror=format-security’ ‘-Wformat’ ‘-Wp,-D_FORTIFY_SOURCE=3’
    ‘-Wp,-D_GLIBCXX_ASSERTIONS’ ‘-march=x86-64’
    ‘-mno-omit-leaf-frame-pointer’

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

# Version 2.1.1
## R CMD check results

0 errors | 0 warnings | 1 note

* This package was archived on 2023-10-13 as it requires the _then_ archived package 'BFpack'.
  This error has now been fixed. 
