# gmx_stacking_NG
NG stands for __Next Generation__

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

This program is designed to work with Gromacs 2018.X versions, 
and we'd highly recommend you stick with this version.

Tests show that it doesn't work with Gromacs 5.X or Gromacs 2019, 
however its compatibility with Gromacs 2016, 2019 is not tested.
