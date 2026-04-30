"""
General utility functions for the project.
"""
import sys
from pathlib import Path


class Tee:
    """Write to both terminal and file simultaneously."""
    def __init__(self, filepath):
        self.file = open(filepath, "w")
        self.terminal = sys.__stdout__

    def write(self, msg):
        self.terminal.write(msg)
        self.file.write(msg)

    def flush(self):
        self.terminal.flush()
        self.file.flush()

    def close(self):
        self.file.close()


def setup_output(filename, output_dir="data/results"):
    """
    Set up dual output to terminal and file.
    Returns the output filepath for reference.

    Usage:
        from src.utils import setup_output, teardown_output
        output_path = setup_output("exploration_output.txt")
        # ... your script ...
        teardown_output(output_path)
    """
    output_path = Path(output_dir) / filename
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sys.stdout = Tee(output_path)
    return output_path


def teardown_output(output_path):
    """
    Restore stdout and print confirmation.
    """
    sys.stdout.close()
    sys.stdout = sys.__stdout__
    print(f"Done — output saved to {output_path}")