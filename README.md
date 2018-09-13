<img src=https://github.com/yufengliang/mbxaspy/blob/master/doc/logo.png>

# Many-Body X-ray Absorption Spectrosopy with PYthon

MBXASPY is a python software package initiated by [Yufeng Liang](https://scholar.google.com/citations?user=xiRU9IEAAAAJ&hl=en) at the [Molecular Foundry (TMF)](http://foundry.lbl.gov), [Lawrence Berkeley National Laboratory (LBNL)](https://www.lbl.gov), for predicting x-ray spectra using the determinant formalism. The determinant formalism is based on the independent-electron approximation as used in the density-functional theory (DFT) and hence the many-body wavefunctions used in Fermi's Golden rule take a form of a single Slater determinant. The orbitals for constructing the initial/final-state Slater determinant are obtained from DFT calculations. The determinant approach for free-electron systems is equivalent to the [MND theory](https://journals.aps.org/rmp/pdf/10.1103/RevModPhys.62.929).

To run MBXASPY, you need to first generate eigenvalues and eigenvectors of your Hamiltonian, either from DFT or from tight-binding models. At this stage, MBXASPY is seamingless interfaced with a Fortran software package ShirleyXAS that performs the DFT calculation at the one-body level. ShirleyXAS is modified based on [Quantum Espresso](https://www.quantum-espresso.org) 4.3 by [David Prendergast](https://scholar.google.com/citations?user=Saf7NMcAAAAJ&hl=en) at TMF. If you would like to perform the ShirleyXAS + MBXASPY calculation for interpreting experimental x-ray spectra from first-principles, please contact [David's group at TMF](http://nanotheory.lbl.gov/people/prendergast.html) for accessiblity and guidance. A manual for performing ShirleyXAS + MBXASPY calculation using an the MBXAS module on the [local cluster of LBNL](http://scs.lbl.gov) can be found in the [doc directory](https://github.com/yufengliang/mbxaspy/blob/master/doc/Manual%20for%20ShirleyXAS%20%2B%20MBXASPY.pdf). Please generously cite the below references if you are using the MBXASPY code.

Currently, you may use MBXASPY to obtain x-ray absorption spectra (XAS) and x-ray photoemission spectra (XPS) at the level of MND theory. One important application so far is to interpret the O *K* edge spectra for transition metal oxides. The authors are also actively looking for other applications that require the determinant method.

The [develop branch](https://github.com/yufengliang/mbxaspy/tree/develop) can handle eigenvalues and eigenvectors from tight-binding models by setting scf_type ='model'. The [RIXS branch](https://github.com/yufengliang/mbxaspy/tree/rixs) for simulating resonant inelastic x-ray scattering (RIXS) and the [emission branch](https://github.com/yufengliang/mbxaspy/tree/emission) for simulaintg x-ray emission spectra based on the determinant method are being developed and tested.

Development
---------
- Interface MBXASPY with other popular DFT codes.

Reference
---------
- Yufeng Liang, John Vinson, Sri Pemmaraju, Walter S. Drisdell, Eric L. Shirley, and David Prendergast,  
*Accurate x-ray spectral predictions: an advanced self-consistent-field approach inspired by many-body perturbation theory*,
[Phys. Rev. Lett. 118, 096402 (2017)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.118.096402).
- Yufeng Liang and David Prendergast,
*Quantum many-body effects in x-ray spectra efficiently computed using a basic graph algorithm*,
[Phys. Rev. B 97, 205127 (2018)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.205127). 
