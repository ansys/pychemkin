PSR network
===========

**PyChemkin** can be used to create an *Equivalence Reactor Network (ERN)* as a reduced-order
(in geometry) model for more complex combustion application such as kiln or combustor of a gas turbine.
There are several ways to build an *ERN*. The first method is to create and solve the reactors one-by-one
starting from the upstream reactor. Adjust/mannipulate the external and outlet streams with the *Mixture*/*Stream*
utilities and apply the desired stream as the inlet for the downstream reactors. The second method is to take advantage
of the *hybridreactornetwork* "model" of **PyChemkin**. Simply defined the reactors with the associated external inlets and
add the reactors to the *hybridreactornetwork*. If there are *recycling* streams from the downstream reactor to the upstream ones,
the *hybridreactornetwork* will solve the entire *ERN* iteratively. In this case, a *tear point* must be explicitly
specified. The last methods is to create the network and its connectivity first, and the entire *ERN* will be solved in a coupled
manner. This methods is the default method in Chemkin GUI, but it is currently under development in **PyChemkin**. Note that currently,
*PSR* is the only reactor model that is allowed in the *ERN*. 

The *PSRChain_xxx* examples show the two methods to build and run a series of linked PSRs (no stream recycling) can be modeled in **PyChemkin**.
The *PSRnetwork* example goes over the steps of creating and running an *ERN* with recycling streams by using the *hybridreactornetwork* method.
