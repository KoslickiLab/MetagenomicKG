import neo4j
import os

# Load environment variables
neo4j_bolt = os.environ['neo4j_bolt']
neo4j_username = os.environ['neo4j_username']
neo4j_password = os.environ['neo4j_password']

# Initialize the Neo4j driver
driver = neo4j.GraphDatabase.driver(neo4j_bolt, auth=(neo4j_username, neo4j_password))

def run_query(query):
    """
    Execute a Cypher query and return the result.
    :param query: a Cypher statement as a string to run
    """
    with driver.session() as session:
        result = session.run(query)
        return list(result)

def node_labels():
    """
    Retrieve distinct node labels from the database.
    """
    labels_query = "MATCH (n) RETURN distinct labels(n) as labels"
    results = run_query(labels_query)
    label_list = [result['labels'][0] for result in results if result['labels']]
    return label_list

def create_index(label_list, property_name):
    """
    Create an index on a given property for each label in the list.
    :param label_list: a list of the node labels in Neo4j
    :param property_name: the property name to index
    """
    for label in label_list:
        # Adjusted to use the new CREATE INDEX syntax
        index_query = f"CREATE INDEX FOR (n:`{label}`) ON (n.{property_name})"
        run_query(index_query)

def constraint(label_list):
    """
    Create a unique constraint on the 'id' property for nodes with the label 'Base'.
    :param label_list: a list of the node labels in Neo4j
    """
    if 'Base' in label_list:
        constraint_query = "CREATE CONSTRAINT ON (n:Base) ASSERT n.id IS UNIQUE"
        run_query(constraint_query)

if __name__ == '__main__':
    node_label_list = node_labels()
    
    # Create Indexes on Node Properties
    create_index(node_label_list, 'node_id')
    create_index(node_label_list, 'node_type')
    create_index(node_label_list, 'all_names')
    create_index(node_label_list, 'description')
    create_index(node_label_list, 'knowledge_source')
    create_index(node_label_list, 'link')
    create_index(node_label_list, 'synonyms')
    create_index(node_label_list, 'is_pathogen')