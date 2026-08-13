SYNTHETIC DATA - DO NOT TRUST ANYTHING IN THIS FOLDER
====================================================

Written by bayesopt_boomerang.m running with
    cfg.solverBackend = 'surrogate'
on 12-Aug-2026 23:33:58.

No COMSOL solve was performed. Every band structure, gap,
fitness value, log row and figure in this folder was produced
by an analytic placeholder (surrogateBoomerangBands), or is a
deliberate failure (the 'stub' backend). The numbers are not
approximations of the real cell - they are not physics at all.

This folder exists only to debug the optimization loop:
objective plumbing, constraints, logging, checkpointing and the
post-run figures, without spending hours of COMSOL time.

Every _bds.mat here additionally carries ds.isSynthetic = true,
and a real run refuses to load a file marked that way. Do not
strip that field, and do not move these files into a real
results folder.
