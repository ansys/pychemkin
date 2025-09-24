Multiprocessing parameter study
===============================

Chemkin application by itself does not support parallel computing, however, in the context of parameter study, each case can run on
its own CPU core to make better use of any available computing power. The examples in this section describe the process of applying Python
multi-threading and multi-processing packages to **PyChemkin** parameter study to improve the overall simulation performance. 
