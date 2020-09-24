# gmx_stacking_fast

## Introduction
This program is a totally rewritten version of the original gmx_stacking program.
As the original program can perform a bit more task than this program 
(non-cg benzene rings with >3 atoms), we choose not to update that program, 
but write a new program (this one) instead.

Purpose of this upgrade is, hopefully, to:
- Earthshaking accelerate the analysis progress
- Enable analysis of very large number of aromatic rings.

These two purposes are very urgent, since the simulation system has been enlarged
greatly in these two years.

## Very important notice