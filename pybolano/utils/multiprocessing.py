import os
from typing import TypedDict

############################################################

__all__ = ["mp_config"]

############################################################


class mp_dict(TypedDict):
    enable: bool
    num_cpus: int
    min_num_args: int

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


############################################################

mp_config = mp_dict()
mp_config["enable"] = True
mp_config["num_cpus"] = os.cpu_count()
mp_config["min_num_args"] = 2
# Skip multiprocessing if the number of elements is small,
# in which case a single core execution is enough.
