# Test exchange proc fucntionality
- checks exchange proc symmetry if one processor has many elements with no MPI sides and other have many MPI sides
  - the one with few MPI sides will not detect all procs but will be detected himself without any problems
- also has mortar elements
- tests many different proc combinations
- Set DO_CORE_SPLIT=ON/OFF to emulate a multi-node environment on a single node
  - either fails on single-node or multi-node (catches different cases)
