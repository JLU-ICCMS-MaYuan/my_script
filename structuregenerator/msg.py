from warnings import warn

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class Formula_Specified_Error(Error):
    """Exception raised for errors in the Compabality.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message