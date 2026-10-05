from importlib.metadata import version

__version__ = version("tarantula_soar")


def hello() -> str:
    return "Hello from tarantula_soar!"
