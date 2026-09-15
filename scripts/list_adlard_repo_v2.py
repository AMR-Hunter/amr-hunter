import json, sys, traceback, urllib.request

API = "https://api.github.com/repos/fowler-lab/tb-bdq-cat"


def get(url):
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "Mozilla/5.0", "Accept": "application/vnd.github+json"},
    )
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))


try:

    def walk(path=""):
        items = get(f"{API}/contents/{path}")
        if isinstance(items, dict):
            print("NOT A LIST:", items, file=sys.stderr)
            return
        for item in items:
            t, name = (item["type"], item["name"])
            if t == "dir":
                print("DIR:", item["path"])
                walk(item["path"])
            else:
                print(
                    f"FILE: {item['path']}  ({item.get('size', '?')} b)  {item.get('download_url')}"
                )

    walk()
except Exception:
    traceback.print_exc()
