# Gaussian_utils
Codes to help dealing with Gaussian Inc. outputs

## hess2freq
it contains and executable that compute the IR frequencies from Hessian Matrix in the file.fchk
To obtain such file, one needs to convert by using **formchk** tool of Gaussian:
```
formchk file.chk file.fchk
```
and then run:
```
hess2freq.exe file.fchk
```

## gau2xyz
Extract optimized geometry from Gaussian Inc. output file (log).
it runs as:
```
gau2xyz.py file.log
```
