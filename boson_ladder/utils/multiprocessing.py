import os

mp_config = {
    "enable" : True,
    "num_cpus" : os.cpu_count(),
    "min_num_args" : 2
    } 
    # Skip multiprocessing if the number of elements is small,
    # in which case a single core execution is enough.