import neo4j
import sys, os
neo4j_bolt = os.getenv('neo4j_bolt')
neo4j_username = os.getenv('neo4j_username')
neo4j_password = os.getenv('neo4j_password')


def run_query(query):
    """
    :param query: a cypher statement as a string to run
    """
    # Start a neo4j session, run a query, then close the session
    session = driver.session()
    query = session.run(query)
    return query


def node_labels():
    # Create a list of dictionaries where each key is "labels(n)"
    # and each value is a list containing a node label
    labels = "MATCH (n) RETURN distinct labels(n)"
    query = run_query(labels)
    label_list = [x.values()[0][0] for x in query]
    return label_list


def create_index(label_list, property_name):
    """
    :param label_list: a list of the node labels in Neo4j
    """
    # For every label in the label list, create an index
    # on the given property name
    for label in label_list:
        index_query = "CREATE INDEX ON :`" + label + "` (" + property_name + ")"
        run_query(index_query)


def constraint(label_list):
    """
    :param label_list: a list of the node labels in Neo4j
    """
    # For every label in the label list, create a unique constraint
    # on the node id property
    constraint_query = "CREATE CONSTRAINT ON (n:Base) ASSERT n.id IS UNIQUE"
    run_query(constraint_query)


if __name__ == '__main__':
    neo4j_password = neo4j_password
    neo4j_user = neo4j_username
    bolt = neo4j_bolt
    driver = neo4j.GraphDatabase.driver(bolt, auth=(neo4j_user, neo4j_password))
    node_label_list = node_labels() + ['Base']

    # Create Indexes on Node Properties
    create_index(node_label_list, "node_type")

    constraint(node_label_list)
    driver.close()