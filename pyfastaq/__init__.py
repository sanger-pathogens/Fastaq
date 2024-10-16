try:
    from importlib.metadata import Distribution
    __version__ = Distribution().from_name('pyfastaq').version
except ImportError:
    from pkg_resources import get_distribution
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
