import json

with open(r"C:\Users\harsh\.gemini\config\plugins\science\skills\literature_search_europepmc\result.json", "r", encoding="utf-8") as f:
    data = json.load(f)

for res in data.get("results", []):
    title = res.get("title", "")
    abstract = res.get("abstractText", "")
    
    if "contact" in title.lower() or "contact" in abstract.lower() or "resistance" in abstract.lower():
        print(f"Title: {title}")
        print(f"Year: {res.get('pubYear')}")
        print(f"Authors: {res.get('authorString')}")
        print(f"DOI: {res.get('doi')}")
        print(f"Abstract snippet: {abstract[:200]}...\n")
