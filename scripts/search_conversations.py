#!/usr/bin/env python3
"""Search the conversation archive bundled with this skill."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
SKILL_DIR = SCRIPT_DIR.parent
INDEX_PATH = SKILL_DIR / "references" / "conversation_index.jsonl"


def normalize(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def snippet(text: str, terms: list[str], limit: int) -> str:
    text = normalize(text)
    lower = text.lower()
    positions = [lower.find(term.lower()) for term in terms if lower.find(term.lower()) >= 0]
    start = max(min(positions) - limit // 3, 0) if positions else 0
    end = min(start + limit, len(text))
    prefix = "..." if start > 0 else ""
    suffix = "..." if end < len(text) else ""
    return prefix + text[start:end].strip() + suffix


def extract_text(content: Any) -> str:
    parts: list[str] = []

    def visit(value: Any) -> None:
        if isinstance(value, str):
            parts.append(value)
        elif isinstance(value, dict):
            for key in ("text", "input_text", "output_text"):
                item = value.get(key)
                if isinstance(item, str):
                    parts.append(item)
            nested = value.get("content")
            if nested is not None:
                visit(nested)
        elif isinstance(value, list):
            for item in value:
                visit(item)

    visit(content)
    return normalize(" ".join(parts))


def load_index(index_path: Path) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    if not index_path.exists():
        return entries
    with index_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                entries.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return entries


def match_score(text: str, terms: list[str]) -> int:
    lower = text.lower()
    score = 0
    for term in terms:
        term_lower = term.lower()
        count = lower.count(term_lower)
        if count == 0:
            return 0
        score += count * max(len(term_lower), 1)
    return score


def index_haystack(entry: dict[str, Any]) -> str:
    pieces: list[str] = [
        str(entry.get("thread_name") or ""),
        str(entry.get("cwd") or ""),
        str(entry.get("source") or ""),
        str(entry.get("originator") or ""),
        str(entry.get("search_text") or ""),
    ]
    for key in ("user_prompts", "assistant_replies"):
        values = entry.get(key) or []
        if isinstance(values, list):
            pieces.extend(str(value) for value in values)
    return "\n".join(pieces)


def search_index(entries: list[dict[str, Any]], terms: list[str], limit: int) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    for entry in entries:
        haystack = index_haystack(entry)
        score = match_score(haystack, terms)
        if not score:
            continue
        results.append(
            {
                "score": score,
                "session_id": entry.get("session_id"),
                "thread_name": entry.get("thread_name"),
                "updated_at": entry.get("updated_at"),
                "cwd": entry.get("cwd"),
                "raw_file": entry.get("raw_file"),
                "snippet": snippet(haystack, terms, 700),
            }
        )
    results.sort(key=lambda item: (item["score"], str(item.get("updated_at") or "")), reverse=True)
    return results[:limit]


def safe_json_lines(path: Path):
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                yield line_number, json.loads(line), line
            except json.JSONDecodeError:
                yield line_number, None, line


def search_deep(
    skill_dir: Path,
    entries: list[dict[str, Any]],
    terms: list[str],
    limit: int,
    context: int,
    include_tools: bool,
) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    seen: set[tuple[str, int]] = set()
    for entry in entries:
        raw_file = entry.get("raw_file")
        if not raw_file:
            continue
        path = skill_dir / str(raw_file)
        if not path.exists():
            continue
        for line_number, obj, raw_line in safe_json_lines(path):
            text = ""
            role = ""
            if isinstance(obj, dict):
                payload = obj.get("payload") or {}
                if obj.get("type") == "response_item" and payload.get("type") == "message":
                    role = str(payload.get("role") or "")
                    text = extract_text(payload.get("content"))
                elif obj.get("type") == "event_msg":
                    text = normalize(str(payload.get("message") or ""))
            if not text and include_tools:
                text = raw_line
            elif not text:
                continue
            score = match_score(text, terms)
            if not score:
                continue
            marker = (str(raw_file), line_number)
            if marker in seen:
                continue
            seen.add(marker)
            results.append(
                {
                    "score": score,
                    "session_id": entry.get("session_id"),
                    "thread_name": entry.get("thread_name"),
                    "updated_at": entry.get("updated_at"),
                    "cwd": entry.get("cwd"),
                    "raw_file": raw_file,
                    "line": line_number,
                    "role": role,
                    "snippet": snippet(text, terms, context),
                }
            )
            if len(results) >= limit:
                return results
    results.sort(key=lambda item: (item["score"], str(item.get("updated_at") or "")), reverse=True)
    return results[:limit]


def print_text(results: list[dict[str, Any]]) -> None:
    if not results:
        print("No matches.")
        return
    for index, item in enumerate(results, start=1):
        title = item.get("thread_name") or "(untitled)"
        print(f"{index}. score={item.get('score')} updated={item.get('updated_at')} title={title}")
        print(f"   session={item.get('session_id')} raw={item.get('raw_file')}")
        if item.get("line"):
            print(f"   line={item.get('line')} role={item.get('role') or ''}")
        if item.get("cwd"):
            print(f"   cwd={item.get('cwd')}")
        print(f"   {item.get('snippet')}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("query", nargs="+", help="One or more search terms. All terms must match.")
    parser.add_argument("--limit", type=int, default=10, help="Maximum results to print.")
    parser.add_argument("--deep", action="store_true", help="Search raw JSONL message lines.")
    parser.add_argument(
        "--include-tools",
        action="store_true",
        help="Include tool calls and outputs in deep search results.",
    )
    parser.add_argument("--context", type=int, default=900, help="Snippet length for deep results.")
    parser.add_argument("--json", action="store_true", help="Print JSONL results.")
    parser.add_argument("--skill-dir", type=Path, default=SKILL_DIR, help="Skill directory.")
    args = parser.parse_args()

    skill_dir = args.skill_dir.resolve()
    entries = load_index(skill_dir / "references" / "conversation_index.jsonl")
    terms = [term for term in args.query if term.strip()]
    if not terms:
        print("No query terms provided.")
        return 2

    if args.deep:
        results = search_deep(
            skill_dir,
            entries,
            terms,
            args.limit,
            args.context,
            args.include_tools,
        )
    else:
        results = search_index(entries, terms, args.limit)

    if args.json:
        for item in results:
            print(json.dumps(item, ensure_ascii=False, sort_keys=True))
    else:
        print_text(results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
