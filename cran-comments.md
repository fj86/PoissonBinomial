## Test environments
* local Manjaro Linux KDE 26.0.2, R 4.5.2
* win-builder (configurations: "release", "oldrelease", "devel")
* mac-builder
* R-hub (configurations: "linux", "m1-san", "macos-arm64", "windows", "valgrind")


## R CMD check results

### local
0 errors | 0 warnings | 0 notes

### win-builder 
0 errors | 0 warnings | 1 note

* checking DESCRIPTION meta-information ... NOTE
  Author field differs from that derived from Authors@R
  
  - must be false alarm, there is no additional Authors field in DESCRIPTION
    file
  - happens only for oldrelease checks, not for release or devel
  
### mac-builder 
0 errors | 0 warnings | 0 notes

### R-hub
0 errors | 0 warnings | 0 notes