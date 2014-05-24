# -*- coding: utf-8 -*-
from __future__ import unicode_literals

"""
Python interface to BioCyc REST API with request caching

Implementation of the BioCyc REST API through the same interface as in Pathway Tools.
Functions are implemented with remote request with cache

Martin Fitzpatrick
May 2014

"""

import os
import errno
import pickle
import requests
import time
from datetime import datetime, timedelta

try:
    import xml.etree.cElementTree as et
except ImportError:
    import xml.etree.ElementTree as et

from .exceptions import BioCycObjectNotFound, BioCycInvalidExpiry, BioCycInvalidDetailLevel

type_converter = {
    'string':str,
    'float':float,
}


from .exceptions import BioCycObjectNotFound
from .singleton import Singleton


DETAIL_NONE = 'none'
DETAIL_LOW = 'low'
DETAIL_FULL = 'full'

DEFAULT_RECORD_EXPIRY = timedelta(weeks=6*4) # Expire after 6 months

DBLINK_URLS = {
    'BIOPATH': "http://www.molecular-networks.com/biopath3/biopath/mols/%s",
    'CAS': "http://www.commonchemistry.org/ChemicalDetail.aspx?ref=%s",
    'CHEBI': "http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:%s",
    'CHEMSPIDER': "http://www.chemspider.com/%s",
    'HMDB': "http://www.hmdb.ca/compounds/%s",
    'KEGG': "http://www.genome.ad.jp/dbget-bin/www_bget?%s",
    'KNAPSACK': "http://kanaya.naist.jp/knapsack_jsp/information.jsp?sname=C_ID&word=%s",
    'LIGAND-CPD': "http://www.genome.ad.jp/dbget-bin/www_bget?%s",
    'NCBI-TAXONOMY-DB': "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=%s",
    'PUBCHEM': "http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=%s",
    'UNIPROT': "http://www.uniprot.org/uniprot/%s",
}

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
        pass


class BioCyc(object):
    __metaclass__ = Singleton

    """
    Basic tools for querying a specific organism via Pathway Tools/BioCyc web API
    """
    _hammer_lock = None
    _hammer_delay = timedelta(seconds=1)

    def __init__(self):
        self.secondary_cache_paths = [] # Not yet implemented
        self.cache_path = os.path.join( os.path.expanduser('~'), '.biocyc' )

        self.set_detail(DETAIL_FULL)
        self.set_organism('HUMAN')
        
        self.expire_records_after = DEFAULT_RECORD_EXPIRY
        
    def set_organism(self, organism):
        self.org_id = organism.upper()
        
    def set_detail(self, detail):
        if detail in [DETAIL_NONE, DETAIL_LOW, DETAIL_FULL]:
            self.detail = detail
        else:
            raise BioCycInvalidDetailLevel
        
    def set_expiry(self, td):
        if type(td) == timedelta:
            self.expire_records_after = td
        else:
            raise BioCycInvalidExpiry
        
    def requestxml(self, url):
    
        # Wait so we don't hammer server
        if self._hammer_lock is not None:
            wait_required = (self._hammer_lock - datetime.now()) + self._hammer_delay
                    
            if not wait_required.days < 0:
                time.sleep(wait_required.seconds)

        self._hammer_lock = datetime.now()
            
        r = requests.get(url)
        if r.status_code == 200:
            # Parse and return the XML
            return et.fromstring(r.text)
        else:
            return False

    def request_api(self, func, org_id, obj):
        return self.requestxml( 'http://websvc.biocyc.org/apixml?fn=%s&id=%s:%s&detail=%s' % (func, org_id, obj, self.detail) )

    def request_obj(self, org_id, obj):
        return self.requestxml( 'http://websvc.biocyc.org/getxml?id=%s:%s&detail=%s' % (org_id, obj, self.detail) )

    def get_from_cache(self, org_id, id):
        '''
        Get an object from the cache
        
        Use all cache folders available (primary first, then secondary in order) and look for the ID in the dir
        if found unpickle and return the object, else return False
        
        FIXME: Check for expiry of object! Return false is expired (will auto-refetch and overwrite)
        '''
        current_time = datetime.now()
        
        for cache in [self.cache_path] + self.secondary_cache_paths:
            read_path = os.path.join( cache, org_id, id )
            try:
                with open(read_path, 'r') as f:
                    obj = pickle.load(f)

            except:
                # Continue to try the next cache
                pass 
                
            else:
                # It worked so we have obj
                # Check for expiry date; if it's not expired return it else continue
                if obj.created_at > current_time - self.expire_records_after:
                    return obj
                    
                # Else continue looking

        # We found nothing (or all expired)
        return None

    def cache(self, obj):
        '''
        Store an object in the cache (this allows temporarily assigning a new cache
        for exploring the DB without affecting the stored version
        '''
        # Check cache path exists for current obj
        write_path = os.path.join( self.cache_path, obj.org_id )
        if not os.path.exists( write_path ): 
            mkdir_p( write_path )

        with open(os.path.join( write_path, obj.id ), 'w') as f:
            pickle.dump( obj, f )
            
    def get(self, ids):
        return self.get_for_org(self.org_id, ids)

    def get_for_org(self, org_id, ids):
        '''
        Returns objects for the given identifiers
        If called with a list returns a list, else returns a single entity
        '''
        t = type(ids)
        if t != list:
            ids = [ids]
            
        objs = []
        
        for id in ids:
            obj = self.get_from_cache(org_id, id)
            if obj is None:
                xml = self.request_obj(org_id, id)
                obj = self.create_obj_from_xml(id, xml)
                self.cache(obj)
            objs.append(obj)
                
        if t != list:
            return objs[0] 
        else:
            return objs
            
            
    def create_obj_from_xml(self, id, xml):
        # Get the object type from the returned XML
        # by matching the provided lists for schema-id
        for o in AVAILABLE_OBJECT_TYPES:
            if o.xml_schema_id:
                x = xml.find(o.xml_schema_id)
                if x: # Found it
                    # Create the base object, populating from the xml
                    # Import the xml to the object (using object specific import_from_xml)
                    obj = o( id=id, from_xml=x)
                    self.cache(obj)
                    return obj
        else:
            raise BioCycObjectNotFound   
            
    def biocyc_obj_url(self, obj):
        return "http://www.biocyc.org/%s/NEW-IMAGE?object=%s"  % (self.org_id, obj)

biocyc = BioCyc()


# Global Pathomx db object class to simplify object display, synonym referencing, etc.
class BioCycEntityBase(object):
    xml_schema_id = None
    ipython_attribs = [
        ('Name', 'name'),
        ('BioCyc ID', 'biocyc_link_html'),
        ('Org ID', 'org_id'),
        ('Synonyms', 'synonyms'),
        ('INCHI', 'inchi'),
        ('Molecular weight', 'molecular_weight'),
        ('Gibbs 0', 'gibbs0'),
        ('Parents', '_parents'),
        ('Instances', '_instances'),
        ('Reactions', '_reactions'),
        ('Pathways', '_pathways'),
        ('Super pathways', '_super_pathways'),
        ('Species', '_species'),
        ('Taxonomic range', '_taxonomic_range'),
        ('Database links', 'dblinks_link_html')]
    
    def __init__(self, id=None, from_xml=None, *args, **kwargs):
        self.type = type(self).__name__.lower()

        self.id = id

        # Parent and child relationships
        self._parents = []
        self._instances = []
        self.name = None
        self.synonyms = []
        
        self.dblinks = {}
        
        # Timestamp object on creation
        self.created_at = datetime.now()
        
        if from_xml:
            self.import_from_xml(from_xml)
        
    def __unicode__(self):
        if self.name:
            return self.name
        else:
            return self.id

    def _str__(self):
        if self.name:
            return self.name
        else:
            return self.id

    def __repr__(self):
        return self.__unicode__()
        
    def _repr_html_(self):
        rows = []
        for l, attr in self.ipython_attribs:
            val = getattr(self, attr, None)
            if val:
                # Manipulate to text
                if type(val) == list:
                    val = ', '.join(val)
                elif type(val) == dict:
                    val = ', '.join( ['%s: %s' % (k,v) for k,v in val.items() ] )
                    
                rows.append('<tr><th>%s</th><td>%s</td></tr>' % (l, val) )
        return "<table>" + ''.join(rows) + "</table>"


    @property
    def cachepath(self):
        return os.path.join(self.org_id, self.type, self.id)

    @property
    def url(self):
        return 'http://biocyc.org/%s/NEW-IMAGE?type=%s&object=%s' % (self.org_id, self.type, self.id)

    def import_from_xml(self, xml):
        '''
        Standard imports for all types of object
        These must fail gracefully, skip if not found
        '''
        self._import_orgid(xml)
        self._import_parents_from_xml(xml)
        self._import_instances_from_xml(xml)
        self._import_common_name(xml)
        self._import_synonyms(xml)
        self._import_dblinks(xml)
        
    def _import_orgid(self, xml):
        self.org_id = xml.attrib['orgid']

    def _import_parents_from_xml(self, xml):
        parents = xml.iterfind('parent')
        for p in parents:
            for o in p:
                # Store a tuple of orgid, identifier
                self._parents.append( o.attrib['frameid'] ) #( o.attrib['orgid'],  ) )

    def _import_instances_from_xml(self, xml):
        instances = xml.iterfind('instance')
        for p in instances:
            for o in p:
                # Store a tuple of orgid, identifier
                self._instances.append( o.attrib['frameid'] ) #( o.attrib['orgid'],  ) )

    def _import_common_name(self, xml):
        e = xml.find('common-name')
        if e is not None:
            self.name = e.text
        
    def _import_synonyms(self, xml):
        es = xml.iterfind('synonym')
        for e in es:
            self.synonyms.append(e.text)
            
        if self.name is None and self.synonyms:
            self.name = self.synonyms[-1] # Apply last synonym if common name not defined
    
    def _import_dblinks(self, xml):
        es = xml.iterfind('dblink')
        for e in es:
            #<dblink-db>LIGAND-CPD</dblink-db><dblink-oid>C00186</dblink-oid>
            self.dblinks[ e.find('dblink-db').text ] = e.find('dblink-oid').text

    def _import_inchi(self, xml):
        self._set_var_from_xml_text( xml, 'inchi', 'inchi') 

    def _import_molecular_weight(self, xml):
        self._set_var_from_xml_text( xml, 'molecular-weight', 'molecular_weight') 

    def _import_gibbs0(self, xml):
        self._set_var_from_xml_text( xml, 'gibbs-0', 'gibbs0' )
    
    def _import_reactions(self, xml):
        self._set_list_ids_from_xml_iter(xml, 'appears-in-right-side-of/Reaction', 'reactions_in_right')
        self._set_list_ids_from_xml_iter(xml, 'appears-in-left-side-of/Reaction', 'reactions_in_left')

    def _import_super_pathways(self, xml):
        self._set_list_ids_from_xml_iter(xml, 'super-pathways/Pathway', '_super_pathways')

    def _import_pathways(self, xml):
        self._set_list_ids_from_xml_iter(xml, 'in-pathway/Pathway', '_pathways')

    def _import_species(self, xml):
        self._set_list_ids_from_xml_iter(xml, 'species/Organism', '_species')
    
    def _import_taxonomic_range(self, xml):
        self._set_list_ids_from_xml_iter(xml, 'taxonomic-range/Organism', '_taxonomic_range')
    
    
    def _set_var_from_xml_text(self, xml, xmlpath, var):
        '''
        Sets a object variable from the xml if it is there
        and passing it through a data conversion based on the variable datatype
        '''
        xmle = xml.find(xmlpath)
        if xmle is not None:
            setattr(self, var, type_converter[ xmle.attrib['datatype'] ]( xmle.text )) 
    
    def _set_list_ids_from_xml_iter(self, xml, xmlpath, var):
        '''
        Set a list variable from the frameids of matching xml entities
        '''
        es = xml.iterfind(xmlpath)
        if es is not None:
            l = []
            for e in es:
                l.append( e.attrib['frameid'] )
            
            setattr(self, var, l)

    @property
    def parents(self):
        return biocyc.get_for_org( self.org_id, self._parents )

    @property
    def instances(self):
        return biocyc.get_for_org( self.org_id, self._instances )

    @property
    def dblinks_link_html(self):
        db = {}
        for k,v in self.dblinks.items():
            if k in DBLINK_URLS:
                db[k] = '<a href="%s">%s</a>' % (DBLINK_URLS[k] % v,v)
            else:
                db[k] = v
        return db

    @property
    def biocyc_link_html(self):
        return '<a href="%s">%s</a>' % (biocyc.biocyc_obj_url(self.id), self.id )



class Compound(BioCycEntityBase):
    xml_schema_id = 'Compound'

    def __init__(self, *args, **kwargs):
        self.inchi = ''
        self.molecular_weight = None
        self.gibbs0 = None
        
        self.reactions_in_right = []
        self.reactions_in_left = []
        
        super(Compound, self).__init__(*args, **kwargs)
    
    def import_from_xml(self, xml):
        super(Compound, self).import_from_xml(xml)
        self._import_inchi(xml)
        self._import_molecular_weight(xml)
        self._import_gibbs0(xml)
        self._import_reactions(xml)
    
                    
    @property
    def _reactions(self):
        return self.reactions_in_right + self.reactions_in_left
       
    @property
    def reactions(self):
        return biocyc.get_for_org( self.org_id, self._reactions )
    
    
class Pathway(BioCycEntityBase):
    xml_schema_id = 'Pathway'

    def __init__(self, *args, **kwargs):
        self._species = []
        self._super_pathways = []
        self._taxonomic_range = []
        super(Pathway, self).__init__(*args, **kwargs)

    def import_from_xml(self, xml):
        super(Pathway, self).import_from_xml(xml)
        self._import_species(xml)
        self._import_super_pathways(xml)
        self._import_taxonomic_range(xml)

    @property
    def species(self):
        return biocyc.get_for_org( self.org_id, self._species )

    @property
    def super_pathways(self):
        return biocyc.get_for_org( self.org_id, self._super_pathways )

    @property
    def taxonomic_range(self):
        return biocyc.get_for_org( self.org_id, self._taxonomic_range )


class Reaction(BioCycEntityBase):
    xml_schema_id = 'Reaction'

    def __init__(self, *args, **kwargs):
        self._pathways = []
        super(Reaction, self).__init__(*args, **kwargs)

    def import_from_xml(self, xml):
        super(Reaction, self).import_from_xml(xml)
        self._import_pathways(xml)

    @property
    def pathways(self):
        return biocyc.get_for_org( self.org_id, self._pathways )
    
class EnzymaticReaction(Reaction):
    xml_schema_id = 'EnzymaticReaction'
    pass

class Protein(BioCycEntityBase):
    xml_schema_id = 'Protein'

class Gene(BioCycEntityBase):
    xml_schema_id = 'Gene'

class DNABindingSite(BioCycEntityBase):
    pass

class Organism(BioCycEntityBase):
    xml_schema_id = 'Organism'
    pass

class Polypeptides(BioCycEntityBase):
    pass

class Promoter(BioCycEntityBase):
    pass

class Complex(BioCycEntityBase):
    pass

class ProteinFeature(BioCycEntityBase):
    pass

class TranscriptionUnit(BioCycEntityBase):
    pass

class tRNA(BioCycEntityBase):
    pass

class Regulation(BioCycEntityBase):
    pass

class Chromosome(BioCycEntityBase):
    pass


AVAILABLE_OBJECT_TYPES = [Compound, Pathway, Reaction, Protein, Gene, DNABindingSite, \
EnzymaticReaction, Organism, Polypeptides, Promoter, Complex, ProteinFeature, \
TranscriptionUnit, tRNA, Regulation]

