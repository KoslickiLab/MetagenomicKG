"""
Utility module for loading configuration from config.yml

This module provides a centralized way to load configuration values
for scripts that need direct access to the configuration file.
"""

import os
import yaml
from typing import Dict, Any

def load_config(config_path: str = None) -> Dict[str, Any]:
    """
    Load configuration from config.yml file.
    
    Args:
        config_path: Path to the config file. If None, looks for config.yml in the root directory.
        
    Returns:
        Dictionary containing the configuration values
        
    Raises:
        FileNotFoundError: If the config file cannot be found
        yaml.YAMLError: If the config file is malformed
    """
    if config_path is None:
        # Look for config.yml in the root directory (parent of build_KG)
        script_dir = os.path.dirname(os.path.abspath(__file__))
        config_path = os.path.join(os.path.dirname(script_dir), 'config.yml')
    
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    try:
        with open(config_path, 'r') as file:
            return yaml.safe_load(file)
    except yaml.YAMLError as e:
        raise yaml.YAMLError(f"Error parsing configuration file: {e}")

def get_build_kg_config(config_path: str = None) -> Dict[str, Any]:
    """
    Get the BUILD_KG_VARIABLES section from the configuration.
    
    Args:
        config_path: Path to the config file. If None, uses default location.
        
    Returns:
        Dictionary containing the BUILD_KG_VARIABLES configuration
    """
    config = load_config(config_path)
    return config.get('BUILD_KG_VARIABLES', {})

def get_neo4j_config(config_path: str = None) -> Dict[str, str]:
    """
    Get Neo4j connection parameters, preferring environment variables over config file.
    
    Args:
        config_path: Path to the config file. If None, uses default location.
        
    Returns:
        Dictionary with neo4j_bolt, neo4j_username, neo4j_password keys
        
    Raises:
        ValueError: If required Neo4j parameters are not found
    """
    build_config = get_build_kg_config(config_path)
    
    neo4j_config = {
        'neo4j_bolt': os.environ.get('neo4j_bolt') or build_config.get('NEO4J_BOLT'),
        'neo4j_username': os.environ.get('neo4j_username') or build_config.get('NEO4J_USERNAME'),
        'neo4j_password': os.environ.get('neo4j_password') or build_config.get('NEO4J_PASSWORD')
    }
    
    missing_params = [key for key, value in neo4j_config.items() if not value]
    if missing_params:
        raise ValueError(f"Missing Neo4j parameters: {', '.join(missing_params)}. "
                        "Set them as environment variables or in config.yml")
    
    return neo4j_config
