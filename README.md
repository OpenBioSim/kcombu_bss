# KCOMBU program

This is a fork of the original KCOMBU program written by

```
==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================
```

[Installation instructions are here](INSTALL.md)

This fork adds code to;

* use `cmake` to support compilation on a range of plaforms
* small fixes to support compilation on, e.g. Windows and MacOS (M1)
* a conda recipe to make it easier to distribute and install binaries

If you use this program please cite the original authors.

1. Kawabata, T. Build-up algorithm for atomic correspondence between chemical structures. J.Chem.Info.Model., 2011,51, 1775-1787.

2. Kawabata T., Nakamura,H. 3D flexible alignment using 2D maximum common substructure: dependence of prediction accuracy on target-reference chemical similarity. J Chem.Info.Model., 2014,54, 1850-1863.

Please also visit the [original website](https://pdbj.org/kcombu) and register
to download the source from them. This way, the original authors can see
who is downloading and using their software.

This fork is not endorsed or supported by the original authors. It has been
made under the terms of the three-clause BSD license specifically to make
the `fkcombu` program easier to install via conda, so that it can be
used as part of [BioSimSpace](https://biosimspace.org).

