lennardjones
============
Simple molecular dynamics simulator, using the Lennard-Jones potential.

This branch uses std::vector<Atom*> to store the atoms in each box, instead of linkedList<Atom*> which the main branch uses. This is mainly for benchmarking and comparison.
