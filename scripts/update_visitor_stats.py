#!/usr/bin/env python3
"""Refresh the country-and-region visitor data used by the homepage map."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from urllib.request import Request, urlopen


COUNTER_ID = "9iO"
COUNTRIES_URL = f"https://s01.flagcounter.com/countries/{COUNTER_ID}/"
METADATA_URL = "https://cdn.jsdelivr.net/npm/world-countries@5.1.0/countries.json"
OUTPUT_PATH = Path(__file__).resolve().parents[1] / "assets/data/visitor-stats.json"

# FlagCounter reports Taiwan separately. The public map presents those visits
# together with China so the ranking, total, fill color, and tooltip stay aligned.
VISITOR_REGION_GROUPS = {
    "CN": ("CN", "China (including Taiwan region)"),
    "TW": ("CN", "China (including Taiwan region)"),
}


class CountryTableParser(HTMLParser):
    """Extract top-level country rows from FlagCounter's public details table."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.rows: list[dict[str, object]] = []
        self._row_stack: list[dict[str, object]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attributes = dict(attrs)
        if tag == "tr":
            self._row_stack.append(
                {"cells": [], "anchors": [], "cell_text": None, "anchor": None}
            )
        elif tag == "td" and self._row_stack:
            self._row_stack[-1]["cell_text"] = []
        elif tag == "a" and self._row_stack:
            self._row_stack[-1]["anchor"] = {
                "href": attributes.get("href") or "",
                "text": "",
            }

    def handle_data(self, data: str) -> None:
        if not self._row_stack:
            return
        row = self._row_stack[-1]
        cell_text = row["cell_text"]
        anchor = row["anchor"]
        if isinstance(cell_text, list):
            cell_text.append(data)
        if isinstance(anchor, dict):
            anchor["text"] += data

    def handle_endtag(self, tag: str) -> None:
        if not self._row_stack:
            return
        row = self._row_stack[-1]
        if tag == "a" and isinstance(row["anchor"], dict):
            anchor = row["anchor"]
            anchor["text"] = " ".join(anchor["text"].split())
            anchors = row["anchors"]
            assert isinstance(anchors, list)
            anchors.append(anchor)
            row["anchor"] = None
        elif tag == "td" and isinstance(row["cell_text"], list):
            cells = row["cells"]
            assert isinstance(cells, list)
            cells.append(" ".join("".join(row["cell_text"]).split()))
            row["cell_text"] = None
        elif tag == "tr":
            completed = self._row_stack.pop()
            completed.pop("cell_text", None)
            completed.pop("anchor", None)
            self.rows.append(completed)


def fetch_text(url: str) -> str:
    request = Request(url, headers={"User-Agent": "Jun-Ma-academic-site/1.0"})
    with urlopen(request, timeout=30) as response:
        return response.read().decode("utf-8")


def parse_country_counts(document: str) -> list[dict[str, object]]:
    parser = CountryTableParser()
    parser.feed(document)
    countries: dict[str, dict[str, object]] = {}

    for row in parser.rows:
        cells = row["cells"]
        anchors = row["anchors"]
        assert isinstance(cells, list)
        assert isinstance(anchors, list)

        country_anchor = next(
            (
                anchor
                for anchor in anchors
                if "/factbook/" in anchor["href"] and f"/{COUNTER_ID}" in anchor["href"]
            ),
            None,
        )
        if country_anchor is None:
            continue

        path_parts = country_anchor["href"].strip("/").split("/")
        if len(path_parts) < 3:
            continue
        source_code = path_parts[1].upper()
        source_name = country_anchor["text"]

        try:
            name_index = cells.index(source_name)
        except ValueError:
            continue

        count = next(
            (int(value) for value in cells[name_index + 1 :] if value.isdigit()),
            None,
        )
        if count is not None:
            code, name = VISITOR_REGION_GROUPS.get(
                source_code, (source_code, source_name)
            )
            if code in countries:
                countries[code]["visitors"] += count
            else:
                countries[code] = {
                    "code": code,
                    "name": name,
                    "visitors": count,
                }

    return sorted(countries.values(), key=lambda country: (-country["visitors"], country["name"]))


def main() -> None:
    country_counts = parse_country_counts(fetch_text(COUNTRIES_URL))
    if not country_counts:
        raise RuntimeError("FlagCounter returned no country rows")

    metadata = json.loads(fetch_text(METADATA_URL))
    metadata_by_iso2 = {country["cca2"]: country for country in metadata}

    for country in country_counts:
        country_metadata = metadata_by_iso2[country["code"]]
        country["iso3"] = country_metadata["cca3"]
        country["map_id"] = country_metadata["ccn3"]

    payload = {
        "updated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "total_visitors": sum(country["visitors"] for country in country_counts),
        "countries": country_counts,
    }

    if OUTPUT_PATH.exists():
        existing = json.loads(OUTPUT_PATH.read_text(encoding="utf-8"))
        if existing.get("countries") == payload["countries"]:
            payload["updated_at"] = existing.get("updated_at", payload["updated_at"])

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
