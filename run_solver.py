from scripts.Solve import solve

# path to the grid network RAW file
casename = 'testcases/GS-4_prior_solution.RAW'

# the settings for the solver
settings = {
    "Tolerance": 1E-05,
    "Max Iters": 1000,
    "Limiting":  False #Voltage limit has been implemented. Please use "True" to enable if needed.
}

# run the solver
solve(casename, settings)