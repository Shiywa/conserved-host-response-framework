#!/usr/bin/env python3
"""Run the adjacent conserved-host-response-framework toolkit without installation."""

from pathlib import Path
import sys

toolkit_root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(toolkit_root / "src"))

from conserved_host_response.cli import main  # noqa: E402


if __name__ == "__main__":
    main()

