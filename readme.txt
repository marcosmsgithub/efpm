Here an algorithmic proposal for the envy-free perfect matching problem is presented.

The maximum weight perfect matching is obtained by the Hungarian method whose code is an adaptation of the one presented by James Payor available at https://github.com/jamespayor/weighted-bipartite-perfect-matching, and then the proposed algorithm for determining the optimal envy-free prices is executed and then compared with the Bellman-Ford algorithm, using the same instances.

Instance sizes can be changed in Line 17 of the efpp.cpp file. The lower and upper bounds of valuations can be changed in Lines 18 and 19, respectively.

Once this is done, just compile the efpp.cpp file and run the executable.
