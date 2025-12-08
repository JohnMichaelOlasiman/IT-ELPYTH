import argparse
import json
import os
import re
import time
from collections import Counter, defaultdict

import requests
from bs4 import BeautifulSoup
try:
    from ord_schema.proto import reaction_pb2
    HAS_ORD = True
except Exception:
    HAS_ORD = False
import base64


BASE_URL = "https://open-reaction-database.org"


def _session():
    s = requests.Session()
    s.headers.update({
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/122 Safari/537.36",
        "Accept": "application/json, text/html;q=0.9",
        "Accept-Language": "en-US,en;q=0.9",
        "Referer": BASE_URL,
    })
    return s


AMINE_KEYWORDS = [
    "amine",
    "aniline",
    "nh2",
    "primary amine",
    "secondary amine",
    "tertiary amine",
]

ARYL_HALIDE_KEYWORDS = [
    "aryl halide",
    "chlorobenz",
    "bromobenz",
    "iodobenz",
    "chloro",
    "bromo",
    "iodo",
    "halide",
]

BASE_KEYWORDS = [
    "carbonate",
    "bicarbonate",
    "hydroxide",
    "tert-butoxide",
    "otbu",
    "triethylamine",
    "dipea",
    "hunig",
    "dbu",
    "dbn",
    "k3po4",
    "k2co3",
    "cs2co3",
    "naoh",
    "koh",
    "naome",
    "kom e",
    "lihmds",
    "nahmds",
    "sodium hydride",
]

METAL_KEYWORDS = [
    "palladium",
    "pd",
    "nickel",
    "ni",
    "copper",
    "cu",
    "iron",
    "fe",
    "ruthenium",
    "ru",
    "iridium",
    "ir",
]

LIGAND_KEYWORDS = [
    "phos",
    "xphos",
    "sphos",
    "johnphos",
    "binap",
    "dppe",
    "dppf",
    "dppp",
    "bipy",
    "phen",
    "bipyridine",
    "phosphine",
    "pcy3",
]

INPUT_KEY_CATEGORY_MAP = {
    "base": "Base",
    "solvent": "Solvent",
    "ligand": "Ligand",
    "metal": "Metal",
    "catalyst": "Metal",
    "amine": "amine",
    "aryl halide": "aryl halide",
    "carboxylic acid": "carboxylic acid",
    "activation agent": "activation agent",
    "additive": "additive",
}


def list_datasets(s: requests.Session):
    u = f"{BASE_URL}/api/datasets"
    r = s.get(u, timeout=30)
    r.raise_for_status()
    return r.json()


def submit_query_for_dataset(s: requests.Session, dataset_id: str, limit: int):
    u = f"{BASE_URL}/api/submit_query"
    params = {"dataset_id": dataset_id, "limit": str(limit)}
    r = s.get(u, params=params, timeout=30)
    r.raise_for_status()
    tid = r.text.strip().strip('"')
    return tid


def fetch_query_result(s: requests.Session, task_id: str, timeout_seconds: int = 120):
    u = f"{BASE_URL}/api/fetch_query_result"
    t0 = time.time()
    delay = 0.5
    while True:
        r = s.get(u, params={"task_id": task_id}, timeout=30)
        if r.status_code == 200:
            return r.json()
        if r.status_code == 202:
            time.sleep(delay)
            delay = min(5.0, delay * 1.5)
        elif r.status_code == 404:
            raise RuntimeError(f"Task {task_id} not found")
        else:
            raise RuntimeError(f"Unexpected status {r.status_code} for task {task_id}")
        if time.time() - t0 > timeout_seconds:
            raise TimeoutError(f"Timeout waiting for task {task_id}")


def get_reaction_summary_html(s: requests.Session, reaction_id: str, compact: bool = False):
    u = f"{BASE_URL}/api/reaction_summary"
    params = {"reaction_id": reaction_id, "compact": "true" if compact else "false"}
    r = s.get(u, params=params, timeout=30)
    r.raise_for_status()
    return r.text


def _pick_name_from_cells(cells):
    for v in cells:
        if not v:
            continue
        if re.fullmatch(r"[0-9 .%]+", v):
            continue
        if v.lower() in ("reagent", "reactant", "solvent", "ligand", "metal", "base"):
            continue
        return v
    return None


def _classify_from_text(name: str, label: str):
    lv = label.lower()
    nm = (name or "").lower()
    if lv == "solvent":
        return "Solvent"
    if lv == "ligand":
        return "Ligand"
    if lv == "metal":
        return "Metal"
    if lv == "base":
        return "Base"
    if lv in ("reagent", "reactant"):
        if any(k in nm for k in BASE_KEYWORDS):
            return "Base"
        if any(k in nm for k in AMINE_KEYWORDS):
            return "amine"
        if any(k in nm for k in LIGAND_KEYWORDS):
            return "Ligand"
        if any(k in nm for k in METAL_KEYWORDS):
            return "Metal"
        if any(k in nm for k in ARYL_HALIDE_KEYWORDS) or ("br" in nm and "c" in nm) or ("cl" in nm and "c" in nm) or ("i" in nm and "c" in nm):
            return "aryl halide"
    return None


def extract_components(html: str):
    soup = BeautifulSoup(html, "lxml")
    out = defaultdict(Counter)
    for tr in soup.find_all("tr"):
        cells = [td.get_text(strip=True) for td in tr.find_all("td")]
        if not cells:
            continue
        label = None
        for v in cells[::-1]:
            lv = v.lower()
            if lv in ("solvent", "reagent", "reactant", "ligand", "metal", "base"):
                label = v
                break
        if not label:
            continue
        name = _pick_name_from_cells(cells[:-1])
        if not name:
            link = tr.find("a")
            name = link.get_text(strip=True) if link else None
        cat = _classify_from_text(name or "", label)
        if cat and name:
            out[cat][name] += 1
    text = soup.get_text(" ", strip=True)
    for lab, cat in [("Base", "Base"), ("Solvent", "Solvent"), ("Ligand", "Ligand"), ("Catalyst", "Metal")]:
        for m in re.finditer(rf"{lab}\s*:\s*([^\n\r<]+)", text, flags=re.IGNORECASE):
            raw = m.group(1)
            parts = [p.strip() for p in re.split(r"[;,]\s*", raw) if p.strip()]
            for p in parts:
                out[cat][p] += 1
    return out


def _enum_name(enum_cls, value: int):
    try:
        return enum_cls.Name(value)
    except Exception:
        return str(value)


def _is_base_ident(value: str):
    v = (value or "").lower()
    if any(k in v for k in BASE_KEYWORDS):
        return True
    if re.search(r"\[o-\].*\[(na|k|li)\+\]", v):
        return True
    if re.search(r"(carbonate|hydroxide|hmds|otbu|tert-?butoxide)", v):
        return True
    return False


def _is_amine_ident(value: str):
    v = (value or "").lower()
    if any(k in v for k in AMINE_KEYWORDS):
        return True
    if re.search(r"n[h]?2", v):
        return True
    return False


def _is_aryl_halide_ident(value: str):
    v = (value or "").lower()
    if any(k in v for k in ARYL_HALIDE_KEYWORDS):
        return True
    if re.search(r"c.*(cl|br|i)", v):
        return True
    return False


def _is_metal_ident(value: str):
    v = (value or "").lower()
    if any(k in v for k in METAL_KEYWORDS):
        return True
    return False


def _is_ligand_ident(value: str):
    v = (value or "").lower()
    if any(k in v for k in LIGAND_KEYWORDS):
        return True
    return False


def aggregate_dataset(s: requests.Session, dataset_id: str, per_dataset_limit: int):
    tid = submit_query_for_dataset(s, dataset_id, per_dataset_limit)
    items = fetch_query_result(s, tid)
    agg_counts = defaultdict(Counter)
    raw = defaultdict(list)
    for item in items:
        rid = item.get("reaction_id")
        if not rid:
            continue
        if HAS_ORD and item.get("proto"):
            try:
                blob = base64.b64decode(item["proto"])
                rxn = reaction_pb2.Reaction()
                rxn.ParseFromString(blob)
                for k, inp in rxn.inputs.items():
                    for comp in inp.components:
                        role_num = getattr(comp, "reaction_role", 0)
                        role_name = _enum_name(reaction_pb2.ReactionRole.ReactionRoleType, role_num) if isinstance(role_num, int) else None
                        idents = []
                        for ident in comp.identifiers:
                            tname = _enum_name(reaction_pb2.CompoundIdentifier.CompoundIdentifierType, ident.type)
                            idents.append({"type": tname, "value": ident.value})
                        cat_by_key = INPUT_KEY_CATEGORY_MAP.get(k.lower())
                        if not cat_by_key and re.fullmatch(r"(?i)m\d+(?:_m\d+)*", k):
                            cat_by_key = k
                        if cat_by_key:
                            for ident in idents:
                                v = ident["value"]
                                raw[cat_by_key].append({"reaction_id": rid, "input_key": k, "reaction_role": role_name, "identifier_type": ident["type"], "value": v})
                                agg_counts[cat_by_key][v] += 1
                        else:
                            for ident in idents:
                                v = ident["value"]
                                if role_name == "SOLVENT":
                                    raw["Solvent"].append({"reaction_id": rid, "input_key": k, "reaction_role": role_name, "identifier_type": ident["type"], "value": v})
                                    agg_counts["Solvent"][v] += 1
                                elif role_name in ("CATALYST",):
                                    if _is_metal_ident(v):
                                        raw["Metal"].append({"reaction_id": rid, "input_key": k, "reaction_role": role_name, "identifier_type": ident["type"], "value": v})
                                        agg_counts["Metal"][v] += 1
                                    if _is_ligand_ident(v):
                                        raw["Ligand"].append({"reaction_id": rid, "input_key": k, "reaction_role": role_name, "identifier_type": ident["type"], "value": v})
                                        agg_counts["Ligand"][v] += 1
                                elif role_name in ("REAGENT", "REACTANT"):
                                    if _is_base_ident(v):
                                        raw["Base"].append({"reaction_id": rid, "input_key": k, "reaction_role": role_name, "identifier_type": ident["type"], "value": v})
                                        agg_counts["Base"][v] += 1
                                    if _is_amine_ident(v):
                                        raw["amine"].append({"reaction_id": rid, "input_key": k, "reaction_role": role_name, "identifier_type": ident["type"], "value": v})
                                        agg_counts["amine"][v] += 1
                                    if _is_aryl_halide_ident(v):
                                        raw["aryl halide"].append({"reaction_id": rid, "input_key": k, "reaction_role": role_name, "identifier_type": ident["type"], "value": v})
                                        agg_counts["aryl halide"][v] += 1
            except Exception:
                pass
        time.sleep(0.1)
    return {"counts": {k: dict(v) for k, v in agg_counts.items()}, "raw": {k: v for k, v in raw.items()}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="*", default=[], help="Dataset IDs to scrape; empty means all")
    ap.add_argument("--limit", type=int, default=50, help="Max reactions per dataset")
    ap.add_argument("--out", type=str, default="ords_components.json", help="Output JSON file")
    ap.add_argument("--timeout", type=int, default=120, help="Timeout per server task")
    ap.add_argument("--base_only", action="store_true")
    ap.add_argument("--smiles_only", action="store_true")
    ap.add_argument("dataset_urls", nargs="*")
    args = ap.parse_args()
    s = _session()
    all_datasets = list_datasets(s)
    ids_from_urls = []
    for u in args.dataset_urls:
        m = re.search(r"ord_dataset-[A-Za-z0-9]+", u)
        if m:
            ids_from_urls.append(m.group(0))
    target_ids = [d["dataset_id"] for d in all_datasets] if not args.datasets and not ids_from_urls else (args.datasets + ids_from_urls)
    result = {}
    for did in target_ids:
        try:
            data = aggregate_dataset(s, did, args.limit)
            if data:
                if args.base_only:
                    base_raw = [x for x in data["raw"].get("Base", []) if (not args.smiles_only or x.get("identifier_type") == "SMILES")]
                    base_counts = data["counts"].get("Base", {})
                    result[did] = {"counts": {"Base": base_counts}, "raw": {"Base": base_raw}}
                else:
                    if args.smiles_only:
                        filtered_raw = {}
                        for k, lst in data["raw"].items():
                            filtered_raw[k] = [x for x in lst if x.get("identifier_type") == "SMILES"]
                        data["raw"] = filtered_raw
                    result[did] = data
        except Exception as e:
            continue
    tmp = json.dumps(result, ensure_ascii=False, indent=2)
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(tmp)
    print(args.out)


if __name__ == "__main__":
    main()

