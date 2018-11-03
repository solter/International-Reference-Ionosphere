# International Reference Ionosphere
This is a fork of the [International Reference Ionosphere (IRI)](http://irimodel.org).
IRI is a emprical model of the ionosphere which provides a reasonable climatological
model of plasma environment in the ionosphere.

The current goals of this repository are to provide the following features:
- [x] Provide an up-to-date git repository that follows the [official release](http://irimodel.org/)
- [ ] Use CI tools to automatically update this to follow the [official release](http://irimodel.org/)
- [ ] Use CI tools to ensure that IRI can be built and start creating tests around it
- [ ] A python interface into IRI to facilitate use in the python
- [ ] Create a set of regression tests to support:
- [ ] A refactorization of IRI into modern fortran, eliminating legacy features such as common blocks and branching gotos

While this ambition is fairly slow going, the master branch
should still serve as a useful git mirror of the [official release](http://irimodel.org/)
in the meantime.
