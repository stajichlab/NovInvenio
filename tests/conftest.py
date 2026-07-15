import sys
from pathlib import Path

# Make the project's `lib/` importable in tests.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'lib'))
