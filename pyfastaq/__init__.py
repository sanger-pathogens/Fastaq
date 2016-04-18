from pkg_resources import get_distribution

try:
    __version__ = get_distribution('pyfastaq').version
except:
    __version__ = 'local'



__all__ = [
    'caf',
    'genetic_codes',
    'utils',
    'sequences',
    'tasks',
    'intervals',
    'runners'
]
from pyfastaq import *
