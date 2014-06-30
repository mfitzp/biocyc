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
import csv
import logging
import re

from datetime import datetime, timedelta
from collections import defaultdict

try:
    import xml.etree.cElementTree as et
except ImportError:
    import xml.etree.ElementTree as et

try:
    # Python 2.x
    import HTMLParser
    html = HTMLParser.HTMLParser()
except:
    # Python 3+
    try:
        import html.parser
        html = html.parser.HTMLParser()
    except:
        # Python 3.5
        import html

strip_tags_re = re.compile(r'<[^>]*?>')
strip_entities_re = re.compile(r'[&;]*') 

type_converter = {
    'string':str,
    'float':float,
    'integer':int,
}


from .exceptions import BioCycObjectNotFound, BioCycInvalidExpiry, BioCycInvalidDetailLevel
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

REACTION_DIRECTIONS = {
    'LEFT-TO-RIGHT': 'forward',
    'RIGHT-TO-LEFT': 'back',
    'REVERSIBLE': 'both',
    'IRREVERSIBLE-LEFT-TO-RIGHT': 'forward',
    'IRREVERSIBLE-RIGHT-TO-LEFT': 'back',
    'PHYSIOL-LEFT-TO-RIGHT': 'forward',
    'PHYSIOL-RIGHT-TO-LEFT': 'back'
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

def to_plain_text(str):
    '''
    Return a plain-text version of a given string

    This is a dumb approach that tags and then removing entity markers
    but this is fine for the content from biocyc where entities are &beta; etc.
    
    Stripping in this way turns these into plaintext 'beta' which is preferable 
    to unicode
    '''
    
    str = strip_tags_re.sub('', str)
    str = strip_entities_re.sub('', str)
    return str

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
        
    def _get_locals(self, table):
        if table in self._locals:
            return self._locals[table]
        else:
            lt = []
            try:
                with open( os.path.join( self.cache_path, self.org_id, table), 'rU') as f:
                    reader = csv.reader(f)
                    for row in reader:
                        lt.append( self.get(row[0]) )
            except:
                return []

            else:
                self._locals[table] = lt
                return lt
                
    @property
    def known_pathways(self):
        return self._get_locals('pathways')

    @property
    def known_genes(self):
        return self._get_locals('genes')

    @property
    def known_compounds(self):
        return self._get_locals('compounds')

    @property
    def known_proteins(self):
        return self._get_locals('proteins')

    @property
    def known_reactions(self):
        return self._get_locals('reactions')
        
    def _get_by_name(self, table, n):
        if table not in self._synonyms:
            nt = {}
            try:
                with open( os.path.join( self.cache_path, self.org_id, table + '-synonyms'), 'rU') as f:
                    reader = csv.reader(f)
                    for row in reader:
                        nt[ row[1] ] = self.get(row[0])
            except Exception as e:
                logging.info(e)
                return None
                
            else:
                self._synonyms[table] = nt


        if n in self._synonyms[table]:
            return self._synonyms[table][n]
        else:
            None

        
    def find_pathway_by_name(self, n):
        return self._get_by_name('pathways', n)

    def find_gene_by_name(self, n):
        return self._get_by_name('genes', n)

    def find_compound_by_name(self, n):
        return self._get_by_name('compounds', n)

    def find_protein_by_name(self, n):
        return self._get_by_name('proteins', n)

    def find_reaction_by_name(self, n):
        return self._get_by_name('reactions', n)

    def find_by_name(self,n ):
        for t in ['pathways', 'genes', 'reactions', 'compounds', 'proteins']:
            o  =  self._get_by_name(t, n)
            if o:
                break
        return o
        
    '''
    This API is incomplete.
    
    The BioCyc remote API for foreignids appears to be broken as of 26.06.2014
    All requests return zero.
    
    def get_via_foreign_id(self, db, id):
        # Implement a cacheing lookup using foreign object IDs
        # Search local table first (CSV-backed) then failing that use a remote
        # lookup. If found, return and get the object, if not return None.
        # Store the lookup for future (?)
        obj = self._get_by_foreign_id(db, id)
        if obj == False: 
            # Not found in cache, and not previously requested
            # None = requested; known to not exist
            # BioCyc uses a plain text API for this bit, random.
            
            r = requests.get('http://websvc.biocyc.org/%s/foreignid?ids=%s' % (db, id))
            if r.status == 200: # OK
                

            wr = [id, obj.id]
            self._foreign_ids[ db ].append( obj )

        
            
            with open( os.path.join( self.cache_path, self.org_id, 'foreign-' + db), 'a+') as f:
                writer = csv.writer(f)
                writer.writerow( wr )

        return obj


    def _get_by_foreign_id(self, db, id):
        if db not in self._foreign_ids:
            nt = {}
            try:
                with open( os.path.join( self.cache_path, self.org_id, 'foreign-' + db), 'rU') as f:
                    reader = csv.reader(f)
                    for row in reader:
                        nt[ row[1] ] = self.get(row[0])
            except Exception as e:
                logging.info(e)
                return None
                
            else:
                self._foreign_ids[db] = nt

        if id in self._foreign_ids[db]:
            return self._foreign_ids[db][id]
        else:
            False
    '''

    def add_to_localstore(self, obj):
        if hasattr(obj, 'localstore'):
            with open( os.path.join( self.cache_path, self.org_id, obj.localstore), 'a+') as f:
                writer = csv.writer(f)
                writer.writerow( [obj.id] )
            self._locals[ obj.localstore ].append( obj )
            
    def add_to_names(self, obj):
        if hasattr(obj, 'localstore'):
        
            name_list = []
            if obj.name is not None:
                name_list.append( obj.name ) # Use plaintext name not the html one
            
            name_list.extend( obj.synonyms )
            name_list = set(name_list) # Only uniques
            
            with open( os.path.join( self.cache_path, self.org_id, obj.localstore + '-synonyms'), 'a+') as f:
                writer = csv.writer(f)
                
                for synonym in name_list:
                    writer.writerow( [obj.id, synonym] )
                    self._synonyms[ obj.localstore ][ synonym ] = obj

    def set_organism(self, organism):
        self.org_id = organism.upper()
        mkdir_p( os.path.join( self.cache_path, self.org_id ) )

        self._locals = defaultdict(list)
        self._synonyms = defaultdict(dict)
        self._foreign_ids = defaultdict(list)
        
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
        
        # Add to localstore (keep track of numbers of objects, etc.)
        self.add_to_localstore(obj)   
        self.add_to_names(obj) 
        
    def get(self, ids, skip_cache=False):
        return self.get_for_org(self.org_id, ids, skip_cache=skip_cache)

    def get_for_org(self, org_id, ids, skip_cache=False):
        '''
        Returns objects for the given identifiers
        If called with a list returns a list, else returns a single entity
        '''
        t = type(ids)
        if t != list:
            ids = [ids]
            
        objs = []
        
        
        for id in ids:
            if id == '' or type(id) is not str: # Empty string
                objs.append(None)
                continue

            if skip_cache ==False:
                obj = self.get_from_cache(org_id, id)
            else:
                obj = None
                
            if obj is None:
                xml = self.request_obj(org_id, id)
                obj = self.create_obj_from_xml(id, xml)
                self.cache(obj) # Will cache either a real object, or a BioCycEntityNotFound
                
            if obj: # Found
                objs.append(obj)
            else:  # Not found (BioCycEntityNotFound)
                objs.append(None)
                
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
                    return obj
        else:
            return BioCycEntityNotFound(id, self.org_id)
            
    def biocyc_obj_url(self, obj):
        return "http://www.biocyc.org/%s/NEW-IMAGE?object=%s"  % (self.org_id, obj)

biocyc = BioCyc()

class BioCycEntityNotFound(object):
    def __init__(self, id=None, org_id=None, *args, **kwargs):
        self.type = type(self).__name__.lower()
        self.id = id
        self.org_id = org_id
        self.created_at = datetime.now()

    def __bool__(self):
        return False
    
    def __nonzero__(self):
        return False

# Global Pathomx db object class to simplify object display, synonym referencing, etc.
class BioCycEntityBase(object):
    xml_schema_id = None
    ipython_attribs = [
        ('Name', 'name_as_html'),
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
        self.name_as_html = None
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

    def __str__(self):
        return self.__unicode__()

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

    def __eq__(self, other):
        return self.type == other.type and self.id == other.id

    def __hash__(self):
        return hash( (self.type, self.id) )

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
            self.name_as_html = e.text
            self.name = to_plain_text(self.name_as_html)
        
    def _import_synonyms(self, xml):
        es = xml.iterfind('synonym')
        for e in es:
            self.synonyms.append(e.text)
            
        if self.name is None and self.synonyms:
            self.name_as_html = self.synonyms[-1] # Apply last synonym if common name not defined
            self.name = to_plain_text(self.name_as_html)
            
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

    def _import_reaction_list(self, xml):
        self._set_list_ids_from_xml_iter(xml, 'reaction-list/Reaction', '_reactions')

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
            setattr(self, var, type_converter[ xmle.attrib.get('datatype', 'string') ]( xmle.text )) 
    
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
            
    def _set_id_from_xml_frameid(self, xml, xmlpath, var):
        '''
        Set a single variable with the frameids of matching entity
        '''
        e = xml.find(xmlpath)
        if e is not None:
            setattr(self, var, e.attrib['frameid'])

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
    localstore = 'compounds'

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

    @property
    def pathways(self):
        return [p for r in self.reactions for p in r.pathways]
    
    
class Pathway(BioCycEntityBase):
    xml_schema_id = 'Pathway'
    localstore = 'pathways'

    def __init__(self, *args, **kwargs):
        self._reactions = []
        self._species = []
        self._super_pathways = []
        self._taxonomic_range = []
        super(Pathway, self).__init__(*args, **kwargs)

    def import_from_xml(self, xml):
        super(Pathway, self).import_from_xml(xml)
        self._import_reaction_list(xml)
        self._import_species(xml)
        self._import_super_pathways(xml)
        self._import_taxonomic_range(xml)

    @property
    def compounds(self):
        return [c for r in self.reactions for c in r.compounds]

    @property
    def reactions(self):
        return biocyc.get_for_org( self.org_id, self._reactions )

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
    localstore = 'reactions'

    def __init__(self, *args, **kwargs):
        self._pathways = []
        self._compounds_left = []
        self._compounds_right = []
        self._enzymatic_reactions = []
        super(Reaction, self).__init__(*args, **kwargs)

    def import_from_xml(self, xml):
        super(Reaction, self).import_from_xml(xml)
        self._import_enzymatic_reaction_objects(xml)
        self._import_pathways(xml)
        self._import_compounds_left(xml)
        self._import_compounds_right(xml)
        self._import_reaction_direction(xml)
        
    @property
    def compounds_left(self):
        return biocyc.get_for_org( self.org_id, self._compounds_left )

    @property
    def compounds_right(self):
        return biocyc.get_for_org( self.org_id, self._compounds_right )

    @property
    def compounds(self):
        return self.compounds_left + self.compounds_right

    @property
    def enzymatic_reactions(self):
        return biocyc.get_for_org( self.org_id, self._enzymatic_reactions )

    @property
    def enzymes(self):
        return [er.enzyme for er in self.enzymatic_reactions]

    @property
    def pathways(self):
        return biocyc.get_for_org( self.org_id, self._pathways )

    def _import_enzymatic_reaction_objects(self, xml):
        # The EnzymaticReaction data in the Reaction XML contains all the information we need
        # to create an EnzymaticReaction object, despite being detail=low.
        # Auto-create them here to avoid unnecessary request
        enzyme_reaction_list = []
        for er in xml.iterfind('enzymatic-reaction/Enzymatic-Reaction'):
            id = er.attrib['frameid']
            enzyme_reaction_list.append(id)
            obj = EnzymaticReaction( id=id, from_xml=er)
            biocyc.cache(obj)
            
        self._enzymatic_reactions = enzyme_reaction_list

    def _import_compounds_left(self, xml):
        self._set_list_ids_from_xml_iter(xml, 'left/Compound', '_compounds_left')

    def _import_compounds_right(self, xml):
        self._set_list_ids_from_xml_iter(xml, 'right/Compound', '_compounds_right')

    def _import_reaction_direction(self, xml):
        self._set_var_from_xml_text( xml, 'reaction-direction', 'direction') 

    
class EnzymaticReaction(BioCycEntityBase):
    xml_schema_id = 'Enzymatic-Reaction'
    localstore = 'enzymaticreactions'

    def __init__(self, *args, **kwargs):
        self._enzyme = None
        self._reaction = None
        super(EnzymaticReaction, self).__init__(*args, **kwargs)

    def import_from_xml(self, xml):
        super(EnzymaticReaction, self).import_from_xml(xml)
        self._import_enzyme(xml)
        self._import_reaction(xml)

    def _import_enzyme(self, xml):
        self._set_id_from_xml_frameid(xml, 'enzyme/Protein', '_enzyme')

    def _import_reaction(self, xml):
        self._set_id_from_xml_frameid(xml, 'reaction/Reaction', '_reaction')

    @property
    def enzyme(self):
        return biocyc.get_for_org( self.org_id, self._enzyme )

    @property
    def reaction(self):
        return biocyc.get_for_org( self.org_id, self._reaction )

    @property
    def pathways(self):
        return self.reaction.pathways


class Protein(BioCycEntityBase):
    xml_schema_id = 'Protein'
    localstore = 'proteins'

    def __init__(self, *args, **kwargs):
        self._gene = None
        self._components = [] # Subunits
        self._catalyzes = []
        self.component_coefficient = None
        super(Protein, self).__init__(*args, **kwargs)

    def import_from_xml(self, xml):
        super(Protein, self).import_from_xml(xml)
        self._import_gene(xml)
        self._import_components(xml)
        self._import_enzymatic_reactions(xml)
        
    @property
    def gene(self):
        return biocyc.get_for_org( self.org_id, self._gene )
    
    def _import_gene(self, xml):
        self._set_id_from_xml_frameid(xml, 'gene/Gene', '_gene')

    @property
    def components(self):
        return biocyc.get_for_org( self.org_id, self._components )

    def _import_components(self, xml):
        self._set_list_ids_from_xml_iter(xml, 'component/Protein', '_components')
        self._set_var_from_xml_text( xml, 'component/coefficient', 'component_coeffecient') 

    @property
    def catalyzes(self):
        return biocyc.get_for_org( self.org_id, self._catalyzes )

    @property
    def reactions(self):
        return [er.reaction for er in self.catalyzes]
    
    @property
    def pathways(self):
        pathway_lists = [er.reaction.pathways for er in self.catalyzes]
        return [p for pl in pathway_lists for p in pl]

    def _import_enzymatic_reactions(self, xml):
        self._set_list_ids_from_xml_iter(xml, 'catalyzes/Enzymatic-Reaction', '_catalyzes')

class Gene(BioCycEntityBase):
    xml_schema_id = 'Gene'
    localstore = 'genes'

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

