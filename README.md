# IRI
This is a fork of the International Reference Ionosphere.
This repository was created to provide the following features:
- [x] Providing an up-to-date git repository that follows the 
[official release](http://irimodel.org/)
- [ ] A refactorization of IRI. This will provide the following advantages:
  * Encapsulation using Fortran 90 modules and derived data types
  * Cleaner code structure through the elimination of labels
    (and thus elimination of gotos)
  * Variable declarations everywhere to avoid implicit types
  * An updated library interface more flexible than the current irisub
- [ ] A python interface into IRI to facilitate use in the python language.

## Download scripts
* download.sh -- Script to download the latest version of iri.
* downloadIndices.sh -- Script to download the lastest indices.
