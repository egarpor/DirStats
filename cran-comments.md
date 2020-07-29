## Test environments

* local R installation, R 3.6.3
* ubuntu 16.04 (on travis-ci), R 3.6.3
* win-builder (release, devel)

Fixed the significant warnings:
    Warning: GNU Extension: Symbol 'nq' is used before it is typed at (1)
    Warning: GNU Extension: Symbol 'kk' is used before it is typed at (1)

One note in Ubuntu Linux 16.04 LTS, R-release, GCC about 
"Found the following (possibly) invalid URLs"
but the cited URLs work properly.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
