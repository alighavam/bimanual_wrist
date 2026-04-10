from __future__ import annotations

import re
import sys
from pathlib import Path


VENDOR_DIR = Path(__file__).resolve().parent / "_vendor"
if VENDOR_DIR.exists():
    sys.path.insert(0, str(VENDOR_DIR))


def _require_reportlab():
    try:
        # type: ignore[import-not-found]
        from reportlab.lib.pagesizes import letter  # noqa: F401
    except Exception as e:  # pragma: no cover
        raise RuntimeError(
            "reportlab is required. Install with:\n"
            "  python -m pip install --target tools/_vendor reportlab"
        ) from e


def _wrap_line(line: str, max_chars: int) -> list[str]:
    if len(line) <= max_chars:
        return [line]
    # Preserve indentation for wrapped continuations.
    indent_match = re.match(r"^(\s+)", line)
    indent = indent_match.group(1) if indent_match else ""
    words = line.strip("\n").split(" ")
    out: list[str] = []
    cur = indent
    for w in words:
        if not w:
            continue
        cand = (cur + (" " if cur.strip() else "") + w) if cur.strip() else (indent + w)
        if len(cand) <= max_chars:
            cur = cand
        else:
            if cur.strip():
                out.append(cur)
            # Start new line; keep indent.
            cur = indent + w
    if cur.strip() or cur == indent:
        out.append(cur)
    return out


def render_text_to_pdf(input_path: Path, output_path: Path) -> None:
    _require_reportlab()
    from reportlab.lib.pagesizes import letter
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.pdfgen import canvas

    # Prefer a mono font for predictable wrapping; fall back to Courier.
    font_name = "Courier"
    try:
        # If a system font is available, user can customize later.
        pass
    except Exception:
        font_name = "Courier"

    page_w, page_h = letter
    margin = 54  # 0.75 inch
    font_size = 10.5
    leading = 14

    c = canvas.Canvas(str(output_path), pagesize=letter)
    c.setTitle(output_path.stem)
    c.setAuthor("Cursor agent")

    c.setFont(font_name, font_size)

    max_width = page_w - 2 * margin

    # Approximate chars-per-line for Courier.
    # Courier is ~0.6 * font_size points per char; this is a practical heuristic.
    approx_char_w = font_size * 0.6
    max_chars = max(40, int(max_width / approx_char_w))

    y = page_h - margin

    text = input_path.read_text(encoding="utf-8", errors="replace").splitlines()
    for raw in text:
        # Turn markdown headings into readable plaintext.
        line = raw.rstrip()
        if line.startswith("## "):
            line = "\n" + line[3:].upper()
        elif line.startswith("# "):
            line = "\n" + line[2:].upper()

        wrapped = _wrap_line(line, max_chars=max_chars) if line else [""]
        for wl in wrapped:
            if y - leading < margin:
                c.showPage()
                c.setFont(font_name, font_size)
                y = page_h - margin
            c.drawString(margin, y, wl)
            y -= leading

    c.save()


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(
            "Usage: python tools/make_pdf_from_text.py <input.md|txt> <output.pdf>",
            file=sys.stderr,
        )
        return 2
    input_path = Path(argv[1]).resolve()
    output_path = Path(argv[2]).resolve()
    if not input_path.exists():
        print(f"Input not found: {input_path}", file=sys.stderr)
        return 2
    output_path.parent.mkdir(parents=True, exist_ok=True)
    render_text_to_pdf(input_path, output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

