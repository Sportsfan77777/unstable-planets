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

    def _check_exception(self, ids, names):
        if len(ids) != len(names):
            raise Exception("The number of names must equal the number of IDs")

    def add(self, id_i, name):
        str_id = str(id_i)

        self._id_dict[str_id] = name
        self._name_dict[name] = id_i

    def add_many(self, ids = [], names = []):
        self._check_exception(ids, names)

        for id_i, name in zip(ids, names):
            self.add(id_i, name)

    def get_id(self, name):
        return self._name_dict[name]

    def get_name(self, id_i):
        str_id = str(id_i)
        return self._id_dict[str_id]

    def read(self, directory = None, name = None):
        if name:
            fn_ids = "ids_%s.p" % name
            fn_names = "names_%s.p" % name
        else:
            fn_ids = "ids.p"
            fn_names = "names.p"

        if directory:
            fn_ids = "%s/%s" % (directory, fn_ids)
            fn_names = "%s/%s" % (directory, fn_names)

        f_ids = open(fn_ids, "rb")
        f_names = open(fn_names, "rb")

        self._id_dict = pickle.load(f_ids)
        self._name_dict = pickle.load(f_names)


    def save(self, directory = None, name = None):
        if name:
            fn_ids = "ids_%s.p" % name
            fn_names = "names_%s.p" % name
        else:
            fn_ids = "ids.p"
            fn_names = "names.p"

        if directory:
            fn_ids = "%s/%s" % (directory, fn_ids)
            fn_names = "%s/%s" % (directory, fn_names)

        f_ids = open(fn_ids, "wb")
        f_names = open(fn_names, "wb")

        pickle.dump(self._id_dict, f_ids)
        pickle.dump(self._name_dict, f_names)

    def find_replace_info_out(self):
        # Writes a new info.out file "named_info.out" using names instead of ids
        f_r = open("info.out", 'r')
        f_w = open("named_info.out", 'w')

        for line in f_r:
            if "ID_" in line:
                id_name = line.split()[0]
                replacement = line.replace(id_name, self.get_name(id_name))
                f_w.write(replacement)
            else:
                f_w.write(line)

        f_r.close()
        f_w.close()

