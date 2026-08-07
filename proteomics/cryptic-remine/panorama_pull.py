#!/usr/bin/env python
"""Minimal PanoramaWeb (LabKey) WebDAV helper: list and download files from a
Panorama folder so a Skyline report can be fed to `skyline-import`.

Panorama exposes every folder over WebDAV at:
    https://panoramaweb.org/_webdav/<Project>/<SubFolder>/@files/

Auth: pass an API key (Panorama: (username) > API Keys > Generate API Key) via
--apikey or the PANORAMA_APIKEY env var. Public folders need no key. The API key
is sent with HTTP Basic auth as user 'apikey'. Nothing is written to disk except
the file you download.

Examples:
    # list a folder
    python panorama_pull.py list  "https://panoramaweb.org/_webdav/MyLab/ALS_CSF/@files/"
    # download one report
    python panorama_pull.py get   "https://panoramaweb.org/_webdav/MyLab/ALS_CSF/@files/report.csv"  ./report.csv

You can also just Map Network Drive to the same _webdav URL in Windows Explorer and
skip this script entirely; this exists for scripted/batch pulls.
"""
import argparse, os, sys
import xml.etree.ElementTree as ET
try:
    import requests
except ImportError:
    sys.exit("pip install requests  (needed for panorama_pull)")


def _auth(apikey):
    apikey = apikey or os.environ.get("PANORAMA_APIKEY")
    return ("apikey", apikey) if apikey else None


def list_folder(url, apikey=None):
    """PROPFIND depth-1 listing; returns list of (href, is_dir, size)."""
    r = requests.request("PROPFIND", url, headers={"Depth": "1"}, auth=_auth(apikey))
    r.raise_for_status()
    ns = {"d": "DAV:"}
    out = []
    for resp in ET.fromstring(r.content).findall("d:response", ns):
        href = resp.findtext("d:href", default="", namespaces=ns)
        is_dir = resp.find(".//d:collection", ns) is not None
        size = resp.findtext(".//d:getcontentlength", default="", namespaces=ns)
        out.append((href, is_dir, size))
    return out


def download(url, dest, apikey=None):
    with requests.get(url, stream=True, auth=_auth(apikey)) as r:
        r.raise_for_status()
        with open(dest, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
    return dest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["list", "get"])
    ap.add_argument("url", help="WebDAV URL (…/@files/… )")
    ap.add_argument("dest", nargs="?", help="local path (for get)")
    ap.add_argument("--apikey", default=None, help="Panorama API key (or PANORAMA_APIKEY env)")
    a = ap.parse_args()
    if a.action == "list":
        for href, is_dir, size in list_folder(a.url, a.apikey):
            print(f"{'d' if is_dir else '-'} {size or '':>12}  {href}")
    else:
        if not a.dest:
            sys.exit("get requires a destination path")
        print("downloaded", download(a.url, a.dest, a.apikey))


if __name__ == "__main__":
    main()
