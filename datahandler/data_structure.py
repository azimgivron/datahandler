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
