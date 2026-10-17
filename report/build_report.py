"""Build a companion document (PDF + Markdown) from a content module in report/ (styles as report/build_report.py).
usage: python report/build_report.py [content_module OUT_STEM]   default: new_biology_content NEW_BIOLOGY
       (run from the repository root)"""
import sys, re, os, importlib
sys.path.insert(0, 'report')
MOD, STEM = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else ("new_biology_content", "report/NEW_BIOLOGY")
C = importlib.import_module(MOD)
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, KeepTogether
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase.pdfmetrics import registerFontFamily
import matplotlib
fd = os.path.join(os.path.dirname(matplotlib.__file__), 'mpl-data', 'fonts', 'ttf')
for nm, f in [('DejaVu', 'DejaVuSans.ttf'), ('DejaVu-Bold', 'DejaVuSans-Bold.ttf'), ('DejaVu-Oblique', 'DejaVuSans-Oblique.ttf'), ('DejaVuMono', 'DejaVuSansMono.ttf')]: pdfmetrics.registerFont(TTFont(nm, os.path.join(fd, f)))
registerFontFamily('DejaVu', normal='DejaVu', bold='DejaVu-Bold', italic='DejaVu-Oblique', boldItalic='DejaVu-Bold')
ss = getSampleStyleSheet()
base = ParagraphStyle('b', parent=ss['BodyText'], fontName='DejaVu', fontSize=9, leading=12.2, spaceAfter=5)
h1 = ParagraphStyle('h1', parent=base, fontName='DejaVu-Bold', fontSize=12.5, leading=15, spaceBefore=10, spaceAfter=6, textColor=colors.HexColor('#1f3b57'))
title = ParagraphStyle('t', parent=base, fontName='DejaVu-Bold', fontSize=17, leading=21, spaceAfter=2)
sub = ParagraphStyle('s', parent=base, fontSize=11, leading=14, textColor=colors.HexColor('#555555'), spaceAfter=10)
cap = ParagraphStyle('c', parent=base, fontSize=7.8, leading=10, textColor=colors.HexColor('#333333'), spaceAfter=10)
cell = ParagraphStyle('cell', parent=base, fontSize=7.4, leading=9.4, spaceAfter=0); cellb = ParagraphStyle('cellb', parent=cell, fontName='DejaVu-Bold')
sumk = ParagraphStyle('sk', parent=base, fontName='DejaVu-Bold', fontSize=8.6, leading=11, textColor=colors.HexColor('#1f3b57'), spaceAfter=0)
sumv = ParagraphStyle('sv', parent=base, fontSize=8.6, leading=11.2, spaceAfter=0)
W = 7.0 * inch
mono = lambda s: s.replace('<font face="Courier">', '<font face="DejaVuMono" size="7">')
def img(path, width=W, max_h=8.4 * inch):
    from PIL import Image as PI
    w, h = PI.open(path).size; ww, hh = width, width * h / w
    if hh > max_h: ww, hh = max_h * w / h, max_h
    return Image(path, width=ww, height=hh)
def pdf(out):
    doc = SimpleDocTemplate(out, pagesize=letter, leftMargin=0.75 * inch, rightMargin=0.75 * inch, topMargin=0.7 * inch, bottomMargin=0.7 * inch, title=C.TITLE, author='')
    st = [Paragraph(C.TITLE, title), Paragraph(C.SUBTITLE + ' · ' + C.DATE, sub)]
    t = Table([[Paragraph(k, sumk), Paragraph(mono(v), sumv)] for k, v in C.SUMMARY], colWidths=[1.35 * inch, W - 1.35 * inch])
    t.setStyle(TableStyle([('BACKGROUND', (0, 0), (-1, -1), colors.HexColor('#f3f6f9')), ('VALIGN', (0, 0), (-1, -1), 'TOP'), ('LINEBELOW', (0, 0), (-1, -2), 0.3, colors.HexColor('#d5dde6')),
                           ('BOX', (0, 0), (-1, -1), 0.6, colors.HexColor('#b8c4d1')), ('TOPPADDING', (0, 0), (-1, -1), 4), ('BOTTOMPADDING', (0, 0), (-1, -1), 4)]))
    st += [Paragraph('<b>Summary</b>', h1), t, Spacer(1, 8)]
    pending = None
    for b in C.BODY:
        if b[0] == 'h1': pending = Paragraph(b[1], h1); continue
        if b[0] == 'fig':
            st.append(KeepTogether(([pending] if pending else []) + [img(b[1]), Paragraph(mono(b[2]), cap)])); pending = None; continue
        if b[0] == 'table':
            data = [[Paragraph(mono(c), cellb if r == 0 else cell) for c in row] for r, row in enumerate(b[1])]
            tt = Table(data, colWidths=[x * inch for x in b[2]], repeatRows=1)
            tt.setStyle(TableStyle([('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#e3e9f0')), ('VALIGN', (0, 0), (-1, -1), 'TOP'), ('GRID', (0, 0), (-1, -1), 0.3, colors.HexColor('#c8d0da')),
                                    ('TOPPADDING', (0, 0), (-1, -1), 3), ('BOTTOMPADDING', (0, 0), (-1, -1), 3)]))
            st.extend([KeepTogether(([pending] if pending else []) + [tt]), Spacer(1, 4)]); pending = None
            continue
        if b[0] == 'p':
            para = Paragraph(mono(b[1]), base)
            st.append(KeepTogether([pending, para]) if pending else para); pending = None
            continue
        if pending: st.append(pending); pending = None
    st += [Paragraph('References', h1)] + [Paragraph(r, base) for r in C.REFS]
    def foot(c, d):
        c.saveState(); c.setFont('DejaVu', 7); c.setFillColor(colors.HexColor('#888888')); c.drawRightString(letter[0] - 0.75 * inch, 0.45 * inch, f'{C.TITLE} · page {d.page}'); c.restoreState()
    doc.build(st, onFirstPage=foot, onLaterPages=foot)
def md(out):
    strip = lambda s: re.sub(r'<font face="Courier">(.*?)</font>', r'`\1`', s).replace('<b>', '**').replace('</b>', '**').replace('<i>', '*').replace('</i>', '*')
    strip2 = lambda s: re.sub(r'<font color="[^"]*">(.*?)</font>', r'\1', strip(s)).replace('&lt;', '<').replace('&gt;', '>')
    L = [f'# {C.TITLE}', f'*{C.SUBTITLE} · {C.DATE}*', '', '## Summary', '']
    for k, v in C.SUMMARY: L.append(f'- **{k}.** {strip2(v)}')
    for b in C.BODY:
        if b[0] == 'h1': L += ['', f'## {b[1]}', '']
        elif b[0] == 'p': L += [strip2(b[1]), '']
        elif b[0] == 'fig': L += [f'![{os.path.basename(b[1])}]({b[3] if len(b) > 3 else b[1]})', '', strip2(b[2]), '']   # optional md-relative path
        elif b[0] == 'table':
            r = b[1]; L.append('| ' + ' | '.join(strip2(c) for c in r[0]) + ' |'); L.append('|' + '---|' * len(r[0]))
            L += ['| ' + ' | '.join(strip2(c).replace('|', '/') for c in x) + ' |' for x in r[1:]]; L.append('')
    L += ['', '## References', ''] + ['- ' + strip2(r) for r in C.REFS]
    open(out, 'w').write('\n'.join(L) + '\n')
if __name__ == '__main__':
    pdf(STEM + '.pdf'); md(STEM + '.md'); print('built', STEM)
