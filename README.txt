The Parallel-Variable Density Approximate Message Passing (P-VDAMP) algorithm 
for reconstructing multi-coil MRI data. This code was used to produce the 
results in 
"Tuning-free multi-coil compressed sensing MRI with Parallel 
Variable Density Approximate Message Passing (P-VDAMP)"
by Charles Millard, Mark Chiew, Jared Tanner, Aaron T Hess and Boris Mailhe. 

*** demos ***

1) PVDAMP_demo.m runs P-VDAMP only. 

2) all_algorithms_demo.m runs P-VDAMP, Optimal FISTA, Shared FISTA, SURE-IT. 
You can also run A-FISTA if you download the authors' code from 
https://gitlab.com/cfmm/datasets/auto-lambda-for-cs-recons 
and add it to the path.  

We have included fully sampled brain k-space data and sensitivity maps for 
the demos. You can choose from acceleration factors R=5,10 for Bernoulli 
or Poisson Disc sampling. The configuration of the algorithm can be 
modified using the opts object.

*** contact *** 

If you have any questions/comments, please feel free to contact Charles 
(Charlie) Millard at <charles.millard@exeter.ox.ac.uk> or Boris Mailhe at
<boris.mailhe@siemens-healthineers.com>

*** copyright and licensing ***

Copyright (C) 2022  Charles Millard
Copyright (C) 2022  Siemens Healthineers

The Software does not have 510K approval,  U.S. Food and Drug Administration (FDA) clearance,
European Medicines Agency (EMA) clearance or CE mark. The Software is intended research purposes only, not for diagnostic or clinical use.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the file GNU_General_Public_License, 
and is also availabe at <https://www.gnu.org/licenses/>
