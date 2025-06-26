PSR network
===========

**PyChemkin** can be used to create an *Equivalence Reactor Network (ERN)* as a reduced-order
(in geometry) model for complex combustion applications such as gas turbine combustor.
Each reactor in the network represents a sub-region (zone) in the actual equipment. Normally these zones
are created according to specific criteria such as temperature, equivalence ratio, or location; and the gas flow
(mean, turbulent, and diffusive) between the zones becomes the internal streams between the reactors.

There are several ways to build an *ERN* in **PyChemkin**. The first method is to create and solve the reactors one-by-one
starting from the most upstream (first) reactor. Adjust/manipulate the external and outlet streams from connected reactors
(that is, solutions from the reactors) using the *Mixture*/*Stream* utilities and apply the desired stream as the inlet for
the downstream reactors. The second method takes advantage of the *hybridreactornetwork* "model" of **PyChemkin**.
Simply defined the reactors with the associated external inlets and add the reactors to the *hybridreactornetwork*.
If there are *recycling* streams from the downstream reactor to the upstream ones, the *hybridreactornetwork* will solve the
entire *ERN* iteratively. In this case, a *tear point* must be explicitly specified. The last method, the coupled method,
is to create the network and its connectivity first, and the entire *ERN* will be solved in a coupled manner. This method
is the default method in Chemkin GUI and is available as the *PSRCluster* "model" in **PyChemkin**.

Chemkin *ERN* has a few limitations:
    * The first reactor of the reactor network must have at least one external inlet.
    * when using the *coupled* model, the entire reactor network can have only one outlet to the surroundings, and the outlet must be attached to
      the last reactor in the network.
    * PSR is the only reactor model allowed in the *PSRCluster* and the *hybridreactornetwork* reactor networks.

The *PSRChain_xxx* examples show the two methods to build and run a series of linked PSRs (no stream recycling) can be modeled in **PyChemkin**.
The *PSRnetwork* example goes over the steps of creating and running an *ERN* with recycling streams by using the *hybridreactornetwork* method.
The *PSRnetwork_coupled* example solves the same *ERN* as the *PSRnetwork* example but utilizes the *coupled* method. Because the sole outlet from
the *coupled ERN* must be attached to the last reactor, the order of the downstream zones is different from the *PSRnetwork* example.
