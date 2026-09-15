import json
import urllib.request

API = "https://api.github.com/repos/fowler-lab/tb-bdq-cat"


def get(url):
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "Mozilla/5.0", "Accept": "application/vnd.github+json"},
    )
    return json.loads(urllib.request.urlopen(req, timeout=60).read().decode("utf-8"))


def walk(path=""):
    for item in get(f"{API}/contents/{path}"):
        t, name = (item["type"], item["name"])
        if t == "dir":
            walk(item["path"])
        else:
            print(
                f"{item['path']}  ({item.get('size', '?')} bytes)  dl={item.get('download_url')}"
            )


walk()
