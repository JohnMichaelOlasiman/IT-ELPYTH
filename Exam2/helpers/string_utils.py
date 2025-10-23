def shout(s: str) -> str:
    """Return s uppercased."""
    if not isinstance(s, str):
        raise TypeError("s must be a string")
    return s.upper()
