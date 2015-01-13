"""
generates an info.out file with different object names than the ones stored by Mercury

Note: This is necessary because Mercury's names are limited to 8 characters. 
More descriptive names are inherently more than 8 characters in most cases.
"""

import sys
from id import ID_Manager

# Initialize ID Manager
id_manager = ID_Manager()
if len(sys.argv) > 1:
	id_manager.read(name = sys.argv[1])
else:
    id_manager.read()

# Run Find+Replace Routine
id_manager.find_replace_info_out()