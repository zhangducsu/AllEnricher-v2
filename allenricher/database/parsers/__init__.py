"""Parsers for supported enrichment database source formats."""

from .go import GOParser
from .kegg import KEGGParser
from .reactome import ReactomeParser
from .do import DOParser
from .disgenet import DisGeNETParser
from .wikipathways import WikiPathwaysParser
from .wikipathways_gpml import WikiPathwaysGPMLParser
from .chea3 import ChEA3Parser
from .trrust import TRRUSTParser
from .htftarget import HTFtargetParser
from .animaltfdb import AnimalTFDBParser

__all__ = ['GOParser', 'KEGGParser', 'ReactomeParser', 'DOParser', 'DisGeNETParser', 'WikiPathwaysParser', 'WikiPathwaysGPMLParser', 'ChEA3Parser', 'TRRUSTParser', 'HTFtargetParser', 'AnimalTFDBParser']
