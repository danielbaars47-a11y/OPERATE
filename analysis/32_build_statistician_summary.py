#!/usr/bin/env python3
"""Build the one-page English statistical handover from reviewed prose.

This is document authoring only. It does not import participant-level files,
reimpute predictors, fit models, or regenerate study results. Use the Codex
bundled Python runtime with python-docx. Render the output separately using
the documents skill render_docx.py and inspect its page PNG before delivery.
"""

from pathlib import Path
import re

from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Mm, Pt, RGBColor


PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "outputs" / "statistician_handover_2026-11-02"
SOURCE = OUT / "OPERATE_methods_and_results.md"
TARGET = OUT / "OPERATE_methods_and_results.docx"


def set_font(style, size, bold=False):
    style.font.name = "Arial"
    style.font.size = Pt(size)
    style.font.bold = bold
    style.font.color.rgb = RGBColor(0, 0, 0)
    fonts = style.element.get_or_add_rPr().get_or_add_rFonts()
    for name in ("ascii", "hAnsi", "eastAsia", "cs"):
        fonts.set(qn(f"w:{name}"), "Arial")


def main():
    text = SOURCE.read_text(encoding="utf-8")
    assert len(text.split()) <= 500, "Condense prose; do not shrink body type."
    document = Document()
    section = document.sections[0]
    section.page_width = Mm(210)
    section.page_height = Mm(297)
    section.top_margin = Mm(18)
    section.bottom_margin = Mm(18)
    section.left_margin = Mm(20)
    section.right_margin = Mm(20)
    section.header_distance = Mm(8)
    section.footer_distance = Mm(8)

    for name in ("Normal", "Title", "Subtitle", "Heading 1", "Heading 2"):
        style = document.styles[name]
        set_font(style, 11)
        style.paragraph_format.space_before = Pt(0)
        style.paragraph_format.space_after = Pt(6)
        style.paragraph_format.line_spacing = 1.04
        style.paragraph_format.widow_control = True
        # Keep titles and headings plain: no inherited theme borders/shading.
        properties = style.element.find(qn("w:pPr"))
        if properties is not None:
            for child in list(properties):
                if child.tag in (qn("w:pBdr"), qn("w:shd")):
                    properties.remove(child)

    set_font(document.styles["Title"], 18, True)
    document.styles["Title"].paragraph_format.space_after = Pt(5)
    set_font(document.styles["Subtitle"], 9)
    document.styles["Subtitle"].paragraph_format.space_after = Pt(11)
    set_font(document.styles["Heading 1"], 11, True)
    document.styles["Heading 1"].paragraph_format.space_before = Pt(5)
    document.styles["Heading 1"].paragraph_format.space_after = Pt(4)
    document.styles["Heading 1"].paragraph_format.keep_with_next = True

    blocks = text.strip().split("\n\n")
    for index, block in enumerate(blocks):
        content = re.sub(r"`([^`]+)`", r"\1", block)
        if content.startswith("# "):
            paragraph = document.add_paragraph(content[2:], style="Title")
            paragraph.paragraph_format.keep_with_next = True
        elif content.startswith("## "):
            document.add_paragraph(content[3:], style="Heading 1")
        elif index == 1:
            document.add_paragraph(content, style="Subtitle")
        else:
            document.add_paragraph(content, style="Normal")

    # English proofing; restrained research memo without a cover or decoration.
    default_rpr = document.styles["Normal"].element.get_or_add_rPr()
    lang = OxmlElement("w:lang")
    lang.set(qn("w:val"), "en-GB")
    default_rpr.append(lang)
    document.core_properties.title = "OPERATE statistical handover"
    document.core_properties.subject = "Existing first-visit QC imputation and prediction results"
    document.core_properties.author = ""
    document.core_properties.last_modified_by = ""
    document.core_properties.comments = ""
    OUT.mkdir(parents=True, exist_ok=True)
    document.save(TARGET)
    print(f"Created {TARGET}; reviewed prose {len(text.split())} words")


if __name__ == "__main__":
    main()
