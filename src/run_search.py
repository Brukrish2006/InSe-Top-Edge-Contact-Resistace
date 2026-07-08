import sys
import subprocess
import os

query = '("InSe" OR "Indium Selenide") AND ("contact resistance" OR "TCAD" OR "FET" OR "transistor") AND (PUB_YEAR:[2020 TO 2026])'

cmd = [
    "uv", "run", "scripts/europepmc_api.py",
    "search",
    query,
    "--max_results", "50",
    "--output", "result.json"
]

subprocess.run(cmd, cwd=r"C:\Users\harsh\.gemini\config\plugins\science\skills\literature_search_europepmc")
