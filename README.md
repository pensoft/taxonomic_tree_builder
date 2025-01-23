# Command Line Argument Parser
This Python script is designed to parse command line arguments for database operations, allowing users to specify table names, database names, and file merging options.

## Features
- Specify the name of the table in the database.
- Set the name of the database with a default option.
- Merge all tables into one.
- Display help information with usage instructions.

## Usage
To use the script, run the following command:

```bash
python arguments.py [OPTIONS]... [FILE]
Options
-t, --table <table_name>: Set the name of the table in the database. Default is auto.
-d, --database <database_name>: Set the name of the database. Default is auto.
-m, --merge: Merge all tables into one.
-h, --help: Display help information with options and arguments.
```

## Example
To run the script with a specific database and table:

```bash
python arguments.py --database my_database --table my_table my_file.txt
```

## Notes
- The filename must be valid; otherwise, a warning will be displayed.
- If the table name is not provided, the user will be prompted to choose from predefined options.

## Merge Command

The `merge` command allows users to combine all specified tables into a single table. This feature is useful for consolidating data from multiple sources or simplifying data management.

### Usage

To use the merge feature, include the `-m` or `--merge` option when running the script. The command will automatically handle the merging process based on the provided database and table parameters.

### Syntax

```bash
python arguments.py --database my_database --table my_table --merge [OPTIONS]
```

## Behavior
- **Merging Process**: When the `--merge` option is specified, the script will:
    - Retrieve all tables from the specified database.
    - Combine the contents of these tables into one unified table.
    - Create a new table with a default name (e.g., `merged_table`), unless specified otherwise.
- **Table Name**: If you want to specify a name for the merged table, use the `--table` option in conjunction with `--merge`.

## Example
To merge all tables and specify a name for the merged table:

```bash
python arguments.py --database my_database --merge --table merged_table
```

## Dependencies
This script requires Python 3.x and the following libraries:

 
`treelib`: Tree structure management.  
`psycopg2`: PostgreSQL database adapter.  
`pick`: Command-line selection tool.  

You can install all requirmenets using pip:

```bash
pip install -r requirements.txt
```
## License
This project is licensed under the MIT License. See the LICENSE file for details.