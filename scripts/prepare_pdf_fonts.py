"""Include every used glyph in embedded Type 1 font metadata."""

import argparse
from pathlib import Path

from pypdf import PdfReader, PdfWriter
from pypdf.generic import (
    ByteStringObject,
    ContentStream,
    NameObject,
    TextStringObject,
)


def prepare_pdf(path):
    reader = PdfReader(path)
    descriptors = {}
    visited = set()

    def scan(stream, resources):
        fonts = resources.get("/Font", {})
        current_font = None
        for operands, operator in ContentStream(stream, reader).operations:
            if operator == b"Tf":
                current_font = fonts.get(operands[0])
            elif operator in (b"Tj", b"TJ", b"'", b'"') and current_font:
                font = current_font.get_object()
                descriptor = font.get("/FontDescriptor")
                encoding = font.get("/Encoding")
                if not descriptor or not encoding:
                    continue
                descriptor = descriptor.get_object()
                encoding = encoding.get_object()
                if font.get("/Subtype") != "/Type1" or "/FontFile" not in descriptor:
                    continue
                if not hasattr(encoding, "get") or "/Differences" not in encoding:
                    continue
                glyphs_by_code = {}
                code = 0
                for item in encoding["/Differences"]:
                    if isinstance(item, int):
                        code = item
                    else:
                        glyphs_by_code[code] = str(item).lstrip("/")
                        code += 1
                key = id(descriptor)
                if key not in descriptors:
                    existing = str(descriptor.get("/CharSet", ""))
                    descriptors[key] = (descriptor, set(existing.split("/")) - {""})
                glyphs = descriptors[key][1]
                texts = operands[0] if operator == b"TJ" else [operands[-1]]
                for text in texts:
                    if isinstance(text, TextStringObject):
                        raw = text.original_bytes
                    elif isinstance(text, ByteStringObject):
                        raw = bytes(text)
                    else:
                        continue
                    glyphs.update(glyphs_by_code[c] for c in raw if c in glyphs_by_code)

        for reference in resources.get("/XObject", {}).values():
            form = reference.get_object()
            if form.get("/Subtype") == "/Form" and id(form) not in visited:
                visited.add(id(form))
                scan(form, form.get("/Resources", resources))

    for page in reader.pages:
        scan(page.get_contents(), page["/Resources"])
    for descriptor, glyphs in descriptors.values():
        descriptor[NameObject("/CharSet")] = TextStringObject(
            "".join("/" + glyph for glyph in sorted(glyphs))
        )
    writer = PdfWriter(clone_from=reader)
    with path.open("wb") as output:
        writer.write(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path, nargs="+")
    for path in parser.parse_args().pdf:
        prepare_pdf(path)
