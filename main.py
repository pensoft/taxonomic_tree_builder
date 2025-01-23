"""
This module processes biological taxonomy data, allowing interaction with a PostgreSQL database.
It reads data from a file, constructs a tree structure of taxonomic data, and facilitates
database operations such as creating tables, inserting data, and indexing.

Dependencies:
- csv
- io
- treelib
- time
- re
- psycopg2
- concurrent.futures
- os
- sys
- getopt
- pick
"""
import csv
import io
from treelib import Tree
import time
import re
import psycopg2 
import concurrent.futures
import os
import sys
import getopt
import pick

start = time.perf_counter()
stime = None
avg_times = []

def key_exists(value, key):
    """
    Check if a key exists in a dictionary.

    Parameters:
    value (dict): The dictionary to check.
    key: The key to check for.

    Returns:
    bool: True if the key exists, False otherwise.
    """
    try:
        value[key]
        return True
    except:
        return False

def get_value_or_default(value, key, default=None):
    """
    Get the value from a dictionary or return a default value.

    Parameters:
    value (dict): The dictionary to check.
    key: The key to look up.
    default: The default value to return if the key is not found.

    Returns:
    The value associated with the key, or the default value.
    """
    if key_exists(value, key):
        return value[key]
    return default

def quote_string_for_db(s):
    """
    Escape a string for safe SQL insertion.

    Parameters:
    s (str): The string to escape.

    Returns:
    str: The escaped string.
    """
    return "\'" + s.replace("'", "''") + "\'"

def is_valid_filename(filename)->bool :
    """
    Check if a filename points to an existing file.

    Parameters:
    filename (str): The path to the file.

    Returns:
    bool: True if the file exists, False otherwise.
    """
    filename = os.path.realpath(filename)
    return os.path.isfile(filename)

def parse_command_line_args():
    """
    Parse command line arguments for the script.

    Returns:
    dict: A dictionary containing database, table name, filename, and merge option.
    """
    def usage(message = None, is_break=False):
        if message != None:
            print("{message}\n".format(message = message))
        if is_break:
            sys.exit(2)
        print("""Usage: python arguments.py [OPTIONS]... [FILE]

Mandatory arguments to long options are mandatory for short options too.
    -h, --host           Set database host.
    -u, --user           Set database user.
    -p, --port           Set database port, default is 5432
    -d, --database       Set name of the database by default is auto.
    -t, --table          Set name of the table in database by default is auto.
    -m, --merge          Merge all tables in one.
    -h, --help           Display help command with options and arguments.
        """)
        sys.exit(2)
        
    try:
        opts, args = getopt.getopt(sys.argv[1:], 't:d:mh', ['table=', 'database=', 'merge', 'help'])
    except:
        usage()

    title = 'Please choose your nomenclature: '
    options = ['NCBI', 'COL', 'GBIF']
    
    database_host = 'taxonomic-db-awt-common.minikube.pensoft.cc'
    database = 'dwca_taxons'
    database_port = 5432
    database_user = 'pensoft'
    
    table = None
    merge = False
    filename = (args[:1] or ['']).pop()
    
    for opt, arg in opts:
        if opt in ('--table', '-t'):
            table = arg
        elif opt in ('--database', '-d'):
            database = arg
        elif opt in ('--port', '-p'):
            database_port = arg
        elif opt in ('--user', '-u'):
            database_user = arg
        elif opt in ('--host', '-h'):
            database_host = arg
        elif opt in ('--merge', '-m'):
            merge = True
        elif opt in ('--help', '-h'):
            usage()
    
    if table is None:
        option, index = pick.pick(options, title, indicator='=>', default_index=1)
        table = "taxon_{0}".format((option or 'default').lower())
    
    if is_valid_filename(filename) == False and merge == False:
        usage('\033[93mWarning: This file is not exists!\033[0m', True)
        
    return {'database': database, 'database_host': database_host, 'database_port': database_port, 'database_user': database_user, 'table': table, 'filename': os.path.realpath(filename), 'merge': merge}
    
class TaxonomicTreeBuilder():
    """
    A class to build and manage a taxonomic tree structure from Checklist Bank DwCA data,
    facilitating interactions with a database.

    Attributes:
        __process_callbacks (list): List of processing callbacks to execute after tree building.
        __extented_callbacks (list): List of extended callbacks for additional processing.
        _connection (object): Database connection object.
        __faild_rows (list): List of rows that failed to process.
        num_read_rows (int): Total number of rows read from the input file.
        rownum (int): Current row number being processed.
        nodes (list): List of node identifiers in the taxonomic tree.
        headers (list): Processed headers for the input data.
        headers_original (list): Original headers from the input data.
        total_bytes (int): Total bytes read from the input file.
        filename (str): Name of the input file.
        table (str): Name of the database table.
    """
    __process_callbacks = []
    __extented_callbacks = []
    _connection = None
    __faild_rows = []
    
    def __init__(self):
        """
        Initializes the TaxonomicTreeBuilder instance,
        creating the tree structure and setting initial attributes.
        """
        self.create_tree()
        self.num_read_rows = 0
        self.rownum = 0
        self.nodes = []
        self.headers = []
        self.headers_original = []
        self.total_bytes = 0
    
    def set_filename(self, filename):
        """Set the filename for input data."""
        self.filename = filename
        
    def get_filename(self):
        """Get the current filename."""
        return self.filename
    
    def set_table(self, table=None):
        """Set the database table name."""
        self.table = table
        
    def get_table(self):
        """
        Get the current table name.
        
        Returns:
        str: The table name.
        """
        if self.table == None:
            table = os.path.splitext(os.path.basename(os.path.realpath(filename)))[0]
            return self.slug(table)
        return self.table
    
    def create_tree(self):
        """Create the initial tree structure with a root node."""
        self.tree = Tree()
        self.tree.create_node('root', 'root')
    
    def add_node_to_tree(self, tag, parent, data):
        """
        Add a node to the tree structure.

        Parameters:
        tag (str): The tag of the node.
        parent (str): The parent node identifier.
        data (dict): Data associated with the node.

        Raises:
        Exception: If there is an error creating the node.
        """
        self.num_read_rows += 1
        identifier = (tag or '').lower()
        parent = (parent or '').lower()
        data.update({'id': self.num_read_rows})
        try:
            self.tree.create_node(tag, identifier, parent, data)
        except Exception as e:
            self.num_read_rows -= 1
            raise Exception(e)
        self.nodes.append(identifier)
    
    def fetch_node(self, nid, extended=False):
        """
        Fetch a node from the tree by its identifier.

        Parameters:
        nid (str): The node identifier.
        extended (bool): Whether to fetch extended data.

        Returns:
        dict: The fetched node or its extended data.
        """
        nid = (nid or '').lower()
        node = self.tree.get_node(nid)
        if extended:
            data = node.data
            for index, row in enumerate(self.read(seek_bytes=data['bytes'])):
                return {'raw':row, 'node': node}
        return node
    
    
    def read(self, header=0, seek_bytes=None, delimiter="\t"):
        """
        Read data from the specified file.

        Parameters:
        header (int): The number of header lines to skip.
        seek_bytes (int): Position to seek in the file.
        delimiter (str): The delimiter used in the file.

        Yields:
        dict: A dictionary containing the row data and metadata.
        """
        with open(self.get_filename(), mode='r') as f:
            line_number = 1
            latest_bytes = None
            
            if seek_bytes is not None:
                f.seek(seek_bytes)
                bytes_position, line = self.read_line(f)
                if not line:  # End of file
                    return
                yield from self.process_line(line, line_number, header, bytes_position, delimiter)
                return

            while True:
                bytes_position, line = self.read_line(f)
                
                if latest_bytes == bytes_position:
                    break
                
                if not line:  # End of file
                    break
                
                yield from self.process_line(line, line_number, header, bytes_position, delimiter)
                latest_bytes = bytes_position
                line_number += 1
    
    def read_line(self, file):
        """Get the current byte position and read a line from the file."""
        bytes_position = file.tell()
        line = file.readline()
        return bytes_position, line
        
    def process_line(self, line, line_number, header, bytes_position, delimiter):
        is_header = header >= line_number
        for row in csv.reader(io.StringIO(line), delimiter=delimiter):
            if is_header:
                self.set_headers(row)
            yield {
                'data': row,
                'bytes': bytes_position,
                'line': line_number,
                'header': is_header
            }
    
    def set_headers(self, row):
        """
        Set the headers for the data.

        Parameters:
        row (list): The row containing headers.
        """
        self.headers_original = row
        self.headers = [self.slug(header) for header in row]
        
    def get_headers(self):
        """Get the processed headers."""
        return self.headers
    
    def combine(self, nid):
        """
        Combine data from a node with its headers.

        Parameters:
        nid (str): The node identifier.

        Returns:
        dict: A dictionary combining headers and data from the node.
        """
        data = self.fetch_node(nid, True)['raw']['data']
        return dict(zip(self.get_headers(), data))
    
    def extends(self, callback: callable):
        """
        Registers a callback for extended processing.

        Args:
            callback (callable): The callback function to register.
        """
        self.__extented_callbacks.append(callback)
        
    def __set_extented_callbacks(self):
        while len(self.__extented_callbacks):
            cb = self.__extented_callbacks.pop()
            self.__process_callbacks.append(cb)
        
    def build_tree(self, row):
        """
        Builds the tree structure from a data row.

        Args:
            row (dict): The data row to process.
        """
        global stime
        if self.rownum % 100001 == 0 or self.rownum == 0:
            stime = time.perf_counter()
        data = row['data']
        taxonId = (data[0] or '')
        parentId = (data[1] or '').lower()
        acceptedNameUsageId = (data[2] or '').lower()
        
        if(row['header']):
            return
        
        if(len(acceptedNameUsageId) > 0):
            try:
                acceptedNameUsageData = self.fetch_node(acceptedNameUsageId, False)
                parentId = self.tree.parent(acceptedNameUsageData.identifier).identifier
            except Exception as e:
                self.__faild_rows.append(row['bytes'])
                return
                
        data_attributes = {
            'bytes': row['bytes']
        }
        if len(parentId):
            try:
                self.add_node_to_tree(taxonId, parentId, data=data_attributes)
            except Exception as e:
                self.__faild_rows.append(row['bytes'])
                return
                
            classification = self.build_classification(taxonId.lower())
            if len(classification) > 0:
                data_attributes.update({'classification': list(dict(classification[:-1]).keys())})
                data_attributes.update({'classification_ids': list(dict(classification[:-1]).values())})
                self.tree.update_node(taxonId.lower(), data=data_attributes)
        else:
            self.add_node_to_tree(taxonId, 'root', data=data_attributes)
        if self.rownum % 100000 == 0 and self.rownum != 0:
            avg_times.append(round(time.perf_counter() - stime, 2))
            print(" Rows:%s, Time:%ss, Avg Time:%ss" % (self.rownum, round(time.perf_counter() - stime, 2), round(sum(avg_times) / len(avg_times), 2)), end="\r")
        self.rownum += 1
    
    def build_classification(self, nid):
        """
        Constructs the classification of a node.

        Args:
            nid (str): The node identifier.

        Returns:
            list: A list of tuples representing the classification.
        """
        def get_parent(nid):
            try:
                node = self.tree.parent(nid)
            except:
                node = None
            return node
        node = get_parent(nid)
        childrens = []
        
        while(node != None):
            nid = node.identifier
            data = node.data
            childrens.append((node.tag, get_value_or_default(data, 'id', 0)))
            node = get_parent(nid)
        return childrens
    
    def processing(self):
        """
        Processes input data to build the taxonomic tree.
        """
        for row in p.read(header=1):
            self.build_tree(row)
        for bytes in self.__faild_rows:
            for row in p.read(header=0, seek_bytes=bytes):
                self.build_tree(row)
        self.__faild_rows = []
        self.__set_extented_callbacks()
        for cb in self.__process_callbacks:
            cb(self, self.nodes)

    def set_connection(self, connection_string):
        """
        Sets the database connection string.

        Args:
            connection_string (str): The connection string for the database.
        """
        self.__connection_string = connection_string
        
    def get_connection_string(self):
        """
        Returns the current database connection string.

        Returns:
            str: The current database connection string.
        """
        return self.__connection_string
            
    def get_connection(self):
        """
        Establishes and returns a database connection.

        Returns:
            object: A connection object to the database.
        """
        return psycopg2.connect(self.__connection_string)

    def create_database(self, database: str = ''):
        """
        Creates a new database with the specified name.

        Args:
            database (str): The name of the database to create.
        """
        connection_string = re.sub(r'(\s?dbname\=\w+\s?)', ' ', self.__connection_string)
        conn = psycopg2.connect(connection_string)
        conn.autocommit = True
        cursor = conn.cursor()
        
        cursor.execute("""
            CREATE DATABASE {} OWNER postgres
        """.format(database))

        cursor.close()
        conn.close()
    
    def get_fields(self) -> dict:
        """
        Returns a dictionary of field names and their data types for table creation.

        Returns:
            dict: A dictionary mapping field names to their data types.
        """
        headers = [x.replace('dwc_', '').replace('col_', '').replace('dcterms_', '') for x in self.get_headers()]
        column_name_with_data_types = dict(map(lambda field: (field, "character varying"), headers))
        column_name_with_data_types.update({"parents": "text[]", "parent_ids": "int[]", "created_at": "timestamp without time zone default NOW()"})
        return dict({**{"id":"bigserial PRIMARY KEY", "label": "character varying"}, **column_name_with_data_types})

    def create_table(self, table_name, column_name_with_data_types = {}):
        """
        Creates a new table in the database.

        Args:
            table_name (str): The name of the table to create.
            column_name_with_data_types (dict): A dictionary of column names and their data types.
        """
        connection = self.get_connection()
        with connection.cursor() as cursor:
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS {0} (
                    {1}
                )
            """.format(table_name, ',\n'.join('\"{}\" {}'.format(key, value) for key, value in column_name_with_data_types.items()) ))
            connection.commit()
        connection.close()
    
    def create_index(self, index, table_name, column_name):
        """
        Creates an index on a specified table column.

        Args:
            index (str): The name of the index to create.
            table_name (str): The name of the table.
            column_name (str): The name of the column to index.
        """
        with self.get_connection() as connection:
            with connection.cursor() as cursor:
                cursor.execute("CREATE INDEX IF NOT EXISTS {0} ON {1}({2})".format(index, table_name, column_name))
    
    def execute_bulk_statements(self, statements, with_records = False):
        """
        Executes multiple SQL statements in bulk.

        Args:
            statements (str): The SQL statements to execute.
            with_records (bool): Whether to return records from the executed statements.

        Returns:
            list: The records returned by the executed statements (if with_records is True).
        """
        connection = self.get_connection()
        connection.autocommit = True
        with connection.cursor() as cursor:
            try:
                cursor.execute(statements)
                if with_records:
                    return cursor.fetchall()
            except Exception as e:
                print(e)
        connection.close()
    
    @staticmethod
    def get_connection_property(property:str, connection_string:str):
        """
        Retrieves a specified property from the connection string.

        Args:
            property (str): The property to retrieve.
            connection_string (str): The connection string to search.

        Returns:
            str: The value of the specified property, or None if not found.
        """
        if property == 'password':
            return
        match = re.search(r'(?<='+property+'\=)(?P<property>\w+)', connection_string)
        if match:
            return match.group('property')

    def slug(self, value):
        """
        Sanitizes a string to create a URL-friendly identifier.

        Args:
            value (str): The string to sanitize.

        Returns:
            str: A sanitized version of the input string.
        """
        return re.sub(r'[\W_]+', '_', str(value)).lower()
        
def save_to_database(app, nodes):
    """
    Saves taxonomic nodes to the database.

    This function creates a database and table if they do not exist, 
    and then inserts the provided taxonomic nodes into the specified table. 
    It utilizes a thread pool for concurrent execution of bulk insert statements to improve performance.
    
    Args:
        app (TaxonomicTreeBuilder): An instance of the TaxonomicTreeBuilder class, 
                                     providing methods for database operations.
        nodes (list): A list of node identifiers to be inserted into the database.
    """
    dbname = app.get_connection_property('dbname', app.get_connection_string())
    try:
        app.create_database(dbname)
    except Exception as e:
        pass
    
    try:
        table_name = app.get_table()
        app.create_table(table_name, app.get_fields())
        app.create_index('taxonid_idx', table_name, 'taxonid')
    except Exception as e:
        print(e)
        
    str_statements = ""
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for i, node in enumerate(nodes):
            node = app.fetch_node(node, True)
            data_attributes = node['node'].data
            data = node['raw']['data']
            id = get_value_or_default(data_attributes, 'id', 1)
            classification = get_value_or_default(data_attributes, 'classification', [])
            classification_ids = get_value_or_default(data_attributes, 'classification_ids', [])
            combined = dict(zip(app.get_headers(), data))
            headers = list(combined.keys())
            values = list(combined.values())
            str_headers = ", ".join(["\"{}\"".format(x.replace('dwc_', '').replace('col_', '').replace('dcterms_', '')) for x in headers + ['parents', 'parent_ids']])
            values_str = ", ".join(["%s" % quote_string_for_db(value) for value in values])
            statement = "INSERT INTO %s (id, %s) VALUES (%d, %s, ARRAY[%s]::text[], ARRAY[%s]::int[]);\n" % (app.table, str_headers, id, values_str, ",".join(["'%s'" % c for c in classification]), ",".join(["%d" % int(c) for c in classification_ids]))
            str_statements += statement
            if i > 0 and i % 100 == 0:
                executor.submit(app.execute_bulk_statements, str_statements)
                str_statements = ""
        if len(str_statements):
            executor.submit(app.execute_bulk_statements, str_statements)
            str_statements = ""

arguments = parse_command_line_args()
database_host = get_value_or_default(arguments, 'database_host')
database = get_value_or_default(arguments, 'database')
database_port = get_value_or_default(arguments, 'database_port')
database_user = get_value_or_default(arguments, 'database_user')
table = get_value_or_default(arguments, 'table')
filename = get_value_or_default(arguments, 'filename')
isMerge = get_value_or_default(arguments, 'merge')

p = TaxonomicTreeBuilder()
p.set_filename(filename)
p.set_table(table)
p.set_connection("host={host} dbname={dbname} user={user} port={port} passfile=./.pgpass".format(host=database_host, dbname=database, user=database_user, port=database_port))

if isMerge:
    p.execute_bulk_statements("""
        DROP TABLE IF EXISTS source_ranking;
        DROP TABLE IF EXISTS taxonranks;
        DROP TABLE IF EXISTS {0};
    """.format(table))
    
    p.create_table('source_ranking', {
        'id': 'bigserial PRIMARY KEY',
        'name': 'character varying',
        'for_zoology': 'integer',
        'for_botany': 'integer',
        'for_mycology': 'integer',
        'general': 'integer',
    })

    p.create_table('taxonranks', {
        'id': 'bigserial PRIMARY KEY',
        'name': 'character varying',
        'ord': 'integer'
    })
    
    p.create_index('idx_taxonranks_name', 'taxonranks', 'name')
    
    p.create_table(table, {
        "id": "bigserial PRIMARY KEY",
        "tid": "integer",
        "taxonid": "character varying",
        "label": "character varying",
        "scientificnameauthorship": "character varying",
        "taxonrank": "character varying",
        "taxonrank_id": "int",
        "taxonomicstatus": "character varying",
        "parents": "text[]",
        "parent_ids": "int[]",
        "source": "character varying",
        "source_id": "integer",
        "kingdom": "character varying"
    })
    
    table_name = p.get_table()
    p.create_index('idx_id', table_name, 'id')
    p.create_index('idx_label', table_name, 'label')
    p.create_index('idx_tid', table_name, 'tid')
    p.create_index('idx_taxonid', table_name, 'taxonid')
    p.create_index('idx_parents', table_name, 'parents')
    p.create_index('idx_taxonrank', table_name, 'taxonrank')
    
    p.execute_bulk_statements("""
        INSERT INTO public.source_ranking (id, name, for_zoology, for_botany, for_mycology, general) VALUES
        (1, 'taxon_worms', 2, 7, 7, 5),
        (8, 'taxon_ipni', 6, 1, 2, 2),
        (3, 'taxon_sfp', 8, 3, 1, 4),
        (6, 'taxon_zoobank', 3, 8, 8, 8),
        (5, 'taxon_ncbi', 5, 6, 6, 7),
        (7, 'taxon_wfo', 7, 2, 3, 3),
        (2, 'taxon_col', 1, 4, 4, 1),
        (4, 'taxon_gbif', 4, 5, 5, 6);
    """)
    
    p.execute_bulk_statements("""
        INSERT INTO public.taxonranks (id, name, ord) VALUES
        (1, 'unranked', 0),
        (2, 'realm', 1),
        (3, 'superkingdom', 5),
        (4, 'kingdom', 10),
        (5, 'subkingdom', 15),
        (6, 'infrakingdom', 20),
        (7, 'division', 25),
        (8, 'subdivision', 30),
        (9, 'superphylum', 35),
        (10, 'infradivision', 35),
        (11, 'parvphylum', 40),
        (12, 'phylum', 45),
        (13, 'subphylum', 50),
        (14, 'infraphylum', 55),
        (15, 'megaclass', 60),
        (16, 'superclass', 65),
        (17, 'gigaclass', 65),
        (18, 'class', 70),
        (19, 'subclass', 75),
        (20, 'subterclass', 80),
        (21, 'infraclass', 85),
        (22, 'parvorder', 90),
        (23, 'superorder', 95),
        (24, 'order', 100),
        (25, 'nanorder', 105),
        (26, 'suborder', 110),
        (27, 'infraorder', 115),
        (28, 'cohort', 120),
        (29, 'subcohort', 125),
        (30, 'superfamily', 130),
        (31, 'epifamily', 135),
        (32, 'family', 140),
        (33, 'subfamily', 145),
        (34, 'infrafamily', 150),
        (35, 'supertribe', 155),
        (36, 'tribe', 160),
        (37, 'subtribe', 165),
        (38, 'infratribe', 170),
        (39, 'suprageneric name', 171),
        (40, 'genus', 175),
        (41, 'subgenus', 180),
        (42, 'supersection botany', 185),
        (43, 'section botany', 190),
        (44, 'subsection botany', 195),
        (45, 'section zoology', 200),
        (46, 'subsection zoology', 205),
        (47, 'superseries', 210),
        (48, 'series', 215),
        (49, 'subseries', 220),
        (50, 'infrageneric name', 225),
        (51, 'species aggregate', 230),
        (52, 'species', 235),
        (53, 'subspecies', 240),
        (54, 'variety', 245),
        (55, 'subvariety', 250),
        (56, 'form', 255),
        (57, 'subform', 260),
        (58, 'forma specialis', 265),
        (59, 'chemoform', 270),
        (60, 'cultivar', 275),
        (61, 'cultivar group', 280),
        (62, 'strain', 285),
        (63, 'morph', 290),
        (64, 'grex', 295),
        (65, 'klepton', 300),
        (66, 'biovar', 305),
        (67, 'pathovar', 310),
        (68, 'chemovar', 315),
        (69, 'natio', 320),
        (70, 'morphovar', 320),
        (71, 'serovar', 321),
        (72, 'proles', 325),
        (73, 'convariety', 330),
        (74, 'mutatio', 335),
        (75, 'lusus', 340),
        (76, 'aberration', 345),
        (77, 'infraspecific name', 350),
        (78, 'infrasubspecific name', 355),
        (79, 'other', 500);
    """)
    
    records = p.execute_bulk_statements("""
        SELECT table_name, '' FROM INFORMATION_SCHEMA.TABLES 
        WHERE table_schema = 'public' AND table_name ILIKE 'taxon\_%'
    """, True)
    
    tables = list(dict(records).keys())
    print("Start processing to merge follow tables:")
    print(",".join(tables))
    for destination_table in tables:
        print("Start to processing table ({0})".format(destination_table))
        print("Start executing update statements in table ({0})".format(destination_table))
        p.execute_bulk_statements("""
            UPDATE {0} SET label = trim(replace(scientificname, scientificnameauthorship, '')); 
        """.format(destination_table))
        connection = p.get_connection()
        connection.autocommit = True
        print("Start executing insert statements in table ({0})".format(table))
        with connection.cursor() as cursor:
            cursor.execute("""
                INSERT INTO {0} 
                (tid, taxonid, label, scientificnameauthorship, taxonrank, taxonomicstatus, parents, parent_ids, source) 
                (SELECT id, taxonid, label, scientificnameauthorship, taxonrank, taxonomicstatus, parents, parent_ids, '{1}' AS source 
                FROM {1})
            """.format(table, destination_table))
        connection.close()
        
    print("Start executing update statements in table ({0})".format(table))
    p.execute_bulk_statements("""
        WITH updated_rows AS (
            SELECT id, name FROM taxonranks
        )
        UPDATE {table}
        SET taxonrank_id = updated_rows.id
        FROM updated_rows
        WHERE {table}.taxonrank = updated_rows.name;
    """.format(table=table))
    
    print("Start executing update statements in table ({0})".format(table))
    p.execute_bulk_statements("""
        WITH updated_rows AS (
            SELECT id, name FROM source_ranking
        )
        UPDATE {table}
        SET source_id = updated_rows.id
        FROM updated_rows
        WHERE {table}.source = updated_rows.name;
    """.format(table=table))
    
    print("".format(table))
    p.execute_bulk_statements("""
        WITH updated_rows AS (
            SELECT 
                ct.id, ct2.label as kingdom 
            FROM cross_taxons as ct
            RIGHT JOIN cross_taxons as ct2 ON 
                ct2.tid = ANY(ct.parent_ids)
                AND ct.source = ct2.source
                AND ct2.taxonrank = 'kingdom'
            ORDER BY
                ct.id ASC
        )
        UPDATE {table}
        SET kingdom = updated_rows.kingdom
        FROM updated_rows
        WHERE {table}.id = updated_rows.id
    """.format(table=table))
    print("Command executed successfully.")
else:
    p.extends(save_to_database)
    p.processing()
    print(round(time.perf_counter() - start, 2), end="\r")
