from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("bootjtk")
except PackageNotFoundError:
    __version__ = "unknown"
