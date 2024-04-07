# Computational Astrophysics N-Body Simulation
This repository contains an implemention for a brute force N-body solver as well as a tree-code implementation made from scratch. The tree-code implementation is based on chapter 4.4 of Völker Springel's "High performance computing and numerical modelling", which can be found in the "files" folder.
A softening factor epsilon for the brute force solver as well as the opening angle for the tree-code implementation can be easily changed in the main file.
The data used contains 50'010 particles that are intended to represent a globular cluster, though the implementation is general enough that it can run for any data that conforms to the same format.


![Teaser1](https://github.com/wwwjones/comp-astro-project/assets/97795524/b985e76f-b462-47bb-aa25-f24126383e8a)


![Teaser2](https://github.com/wwwjones/comp-astro-project/assets/97795524/1894648d-2e45-43eb-a87e-5998f43d1dcc)


![Teaser3](https://github.com/wwwjones/comp-astro-project/assets/97795524/86254fec-0afc-43b1-8e29-03b867001d78)


This project uses Eigen as a dependency and the cmake file assumes that the Eigen installation is included on the PATH.

This project was initially made as the final project of the ETH Zürich AST 245.1 Computational Astrophysics course.
