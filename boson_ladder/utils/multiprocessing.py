import os

class CustomDict(dict):
    def __setitem__(self, key, value):
        valid_keys = ["enable", "num_cpus", "min_num_args"]
        if key not in valid_keys:
            msg = f"The key [{key}] is not valid. Valid keys: {valid_keys}."
            raise KeyError(msg)

        if key == "min_num_args":
            if value >= 2:
                super().__setitem__(key, value)
            else:
                super().__setitem__(key, 2)
        else:
            super().__setitem__(key, value)

mp_config = CustomDict()
mp_config["enable"] = True
mp_config["num_cpus"] = os.cpu_count()
mp_config["min_num_args"] = 2
    # Skip multiprocessing if the number of elements is small,
    # in which case a single core execution is enough.