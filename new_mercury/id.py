"""
ID class to convert stupid mercury names to actual names that aren't limited by 8 characters
"""

import pickle

class ID_Manager():
    """ 
    (1) stores (id, name) pairs retrievable through id or name
    (2) includes read and save options
    """

    def __init__(self, ids = [], names = []):
        self._check_exception(ids, names)

        self._id_dict = {}
        self._name_dict = {}

        self.add_many(ids, names)

    def _check_exception(ids, names):
        if len(ids) != len(names):
            raise Exception("The number of names must equal the number of IDs")

    def add(id_i, name):
        str_id = str(id_i)

        self._id_dict[name] = id_i
        self._name_dict[str_id] = name

    def add_many(ids = [], names = []):
        self._check_exception(ids, names)

        for id_i, name in zip(ids, names):
            self.add(id_i, name)

    def get_id(name):
        return self._id_dict[name]

    def get_name(id_i):
        str_id = str(id_i)
        return self._id_dict[str_id]

    def read(self, name = None):
        if name:
            fn_ids = "ids_%s.p" % name
            fn_names = "names_%s.p" % name
        else:
            fn_ids = "ids.p"
            fn_names = "names.p"

        f_ids = open(fn_ids, "rb")
        f_names = open(fn_names, "rb")

        self._id_dict = pickle.load(f_ids)
        self._name_dict = pickle.load(f_names)


    def save(name = None):
        if name:
            fn_ids = "ids_%s.p" % name
            fn_names = "names_%s.p" % name
        else:
            fn_ids = "ids.p"
            fn_names = "names.p"

        f_ids = open(fn_ids, "wb")
        f_names = open(fn_names, "wb")

        pickle.dump(f_ids, self._id_dict)
        pickle.dump(f_names, self._name_dict)

