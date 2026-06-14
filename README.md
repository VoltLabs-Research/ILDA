# Interfacial Line Defect Analysis (ILDA)

> This is a C++23 port of the original ILDA algorithm.
> The original is an OVITO Pro Python modifier (https://gitlab.com/mmod_public/ilda/-/blob/release/ILDA/utilities.py?ref_type=heads).

## Description

ILDA is designed to have similar functionality to the well-established dislocation extraction algorithm (DXA). Given an atomistic dataset containing a semi-coherent interface, ILDA will identify all interfacial dislocations and disconnections, and insert a line segment representation of the defects. Each line segment has a Burgers vector and step height associated with it. ILDA can be utilized in conjunction with DXA to provide a comprehensive analysis of line defects in an atomistic system.

## Tutorial
First time users may find [this tutorial](https://www.youtube.com/watch?v=O0tElGuW7pA) useful, which demonstrates how to use ILDA from start-to-finish. Note that this video describes how to run ILDA as a Python script.

## Authors and acknowledgment
ILDA was developed by Dr. Nipal Deka, Prof. Ryan Sills, and Dr. Inam Lalani of the Department of Materials Science and Engineering at Rutgers University, in collaboration with Dr. Alexander Stukowski of OVITO GmbH. Its development was supported by the U.S. Department of Energy, Office of Science, Basic Energy Sciences, under Award # DE-SC0022154 (N.D., R.B.S., and I.L.). The ILDA algorithm is documented in this [journal article](https://doi.org/10.1016/j.actamat.2023.119096).

## License
ILDA is open-source software, see the LICENSE file for more details.