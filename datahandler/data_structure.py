import os
from dataclasses import dataclass

@dataclass
class NetworkConfig:
    """
    Generic configuration for downloading a network file.

    Attributes:
        url (str): The download URL.
        filename (str): The local filename to save the retrieved file.
    """

    url: str
    filename: str


@dataclass
class OMIMAPIConfig:
    """
    Configuration for the OMIM REST API.

    Attributes:
        api_key (str): OMIM API key (from OMIM_API_KEY env var).
        base_url (str): Base endpoint for OMIM API calls.
    """

    api_key: str = os.getenv("OMIM_API_KEY")
    base_url: str = "https://api.omim.org/api"