from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

try:
    import yaml as _yaml
except ImportError:  # optional: the bundled fallback supports the toolkit's config subset
    _yaml = None


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)


def read_numeric_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep)


def write_table(frame: pd.DataFrame, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    frame.to_csv(path, sep=sep, index=False)
    return path


def load_yaml(path: str | Path) -> dict[str, Any]:
    text = Path(path).read_text(encoding="utf-8")
    if text.lstrip().startswith(("{", "[")):
        value = json.loads(text)
    elif _yaml is not None:
        value = _yaml.safe_load(text)
    else:
        value = _minimal_yaml_load(text)
    return value or {}


def write_yaml(value: Any, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if _yaml is not None:
        text = _yaml.safe_dump(value, sort_keys=False, allow_unicode=True)
    else:
        text = _minimal_yaml_dump(value)
    path.write_text(text, encoding="utf-8")
    return path


def write_json(value: Any, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False, sort_keys=True)
        handle.write("\n")
    return path


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"1", "true", "yes", "y", "t"})


def _yaml_scalar(value: str):
    value = value.strip()
    if not value:
        return ""
    if (value.startswith('"') and value.endswith('"')) or (value.startswith("'") and value.endswith("'")):
        return value[1:-1]
    lowered = value.lower()
    if lowered in {"true", "false"}:
        return lowered == "true"
    if lowered in {"null", "none", "~"}:
        return None
    if value.startswith("[") and value.endswith("]"):
        body = value[1:-1].strip()
        return [] if not body else [_yaml_scalar(part) for part in body.split(",")]
    if re.fullmatch(r"-?\d+", value):
        return int(value)
    if re.fullmatch(r"-?(?:\d+\.\d*|\d*\.\d+)", value):
        return float(value)
    return value


def _minimal_yaml_load(text: str):
    raw_lines = []
    for raw in text.splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        indent = len(raw) - len(raw.lstrip(" "))
        raw_lines.append((indent, raw.strip()))

    def parse_block(index: int, indent: int):
        is_list = raw_lines[index][1].startswith("- ")
        container = [] if is_list else {}
        while index < len(raw_lines):
            current_indent, content = raw_lines[index]
            if current_indent < indent:
                break
            if current_indent > indent:
                raise ValueError(f"Unsupported YAML indentation near: {content}")
            if is_list:
                if not content.startswith("- "):
                    break
                container.append(_yaml_scalar(content[2:]))
                index += 1
                continue
            if ":" not in content:
                raise ValueError(f"Expected YAML key:value near: {content}")
            key, rest = content.split(":", 1)
            key = key.strip()
            rest = rest.strip()
            if rest:
                container[key] = _yaml_scalar(rest)
                index += 1
            else:
                if index + 1 < len(raw_lines) and raw_lines[index + 1][0] > indent:
                    child_indent = raw_lines[index + 1][0]
                    child, index = parse_block(index + 1, child_indent)
                    container[key] = child
                else:
                    container[key] = {}
                    index += 1
        return container, index

    if not raw_lines:
        return {}
    value, end = parse_block(0, raw_lines[0][0])
    if end != len(raw_lines):
        raise ValueError("Could not parse complete YAML document")
    return value


def _minimal_yaml_dump(value: Any, indent: int = 0) -> str:
    lines = []
    prefix = " " * indent
    if isinstance(value, dict):
        for key, child in value.items():
            if isinstance(child, (dict, list)):
                lines.append(f"{prefix}{key}:")
                lines.extend(_minimal_yaml_dump(child, indent + 2).rstrip().splitlines())
            else:
                lines.append(f"{prefix}{key}: {_dump_scalar(child)}")
    elif isinstance(value, list):
        lines.extend(f"{prefix}- {_dump_scalar(child)}" for child in value)
    else:
        lines.append(f"{prefix}{_dump_scalar(value)}")
    return "\n".join(lines) + "\n"


def _dump_scalar(value: Any) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)
