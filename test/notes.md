# Notes on code optimization

- Initial benchmarks
    - ODE: 6.9ms, approx 1 MB
    - IBM: 20s, 3.7 GB (t_max = 56)
- replacing mapslices with reduce
    - ODE: 4.5 ms, approx 1 MB
    - IBM: 3.7 GB allocations
    - ODE got a little faster, but no big change for IBM
    - wait...profiler still shows mapslices where we should have reduce
    - re-starting session + adjusting
        - ODE: 4.8ms, 1.2 MB
        - IBM: 20s
- removing h_z = sum(@...)
    - IBM: 17-18sm 3.4 GB
    - So some improvement...
- adding for-loop for h_z
    - IBM: 20s, 3.4 GB
- inlining drc functs
    - IBM: 15s, 3.0 GB
- fastmath f_X
    -
- fastmath does not seem to do much good
- TODO: transpose C_W matrices (column major)