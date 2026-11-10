from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from lxml import etree as E
from copy import deepcopy
from openpyxl import load_workbook
P=Path('/tmp/nnrti-nov4-stage/Supplementary-Table-7.xlsx');F=P.parent/'table7-formatted.xlsx';N='http://schemas.openxmlformats.org/spreadsheetml/2006/main';ns={'s':N};REL='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
with ZipFile(P) as z,ZipFile(F) as f:
 updates={};styles=E.fromstring(z.read('xl/styles.xml'));xfs=styles.find('{'+N+'}cellXfs');newstyle=len(xfs);xf=deepcopy(xfs[1]);align=xf.find('{'+N+'}alignment')
 if align is None:align=E.SubElement(xf,'{'+N+'}alignment')
 align.set('wrapText','1');xf.set('applyAlignment','1');xfs.append(xf);xfs.set('count',str(len(xfs)));updates['xl/styles.xml']=E.tostring(styles)
 book=E.fromstring(f.read('xl/workbook.xml'));rels=E.fromstring(f.read('xl/_rels/workbook.xml.rels'));targets={x.get('Id'):x.get('Target') for x in rels}
 for number,name in [(7,'NNRTI tyrosines'),(8,'Tyrosine atoms')]:
  item=next(s for s in book.find('{'+N+'}sheets') if s.get('name')==name);path=targets[item.get('{'+REL+'}id')];path=path.lstrip('/') if path.startswith('/') else 'xl/'+path
  authored=E.fromstring(f.read(path));root=E.fromstring(z.read(f'xl/worksheets/sheet{number}.xml'))
  cols=deepcopy(authored.find('{'+N+'}cols'));root.insert(list(root).index(root.find('{'+N+'}sheetData')),cols)
  first=root.find('{'+N+'}sheetData')[0];first.set('ht','72' if number==7 else '28');first.set('customHeight','1')
  for c in first:c.set('s',str(newstyle))
  updates[f'xl/worksheets/sheet{number}.xml']=E.tostring(root)
 tmp=P.with_suffix('.tmp.xlsx')
 with ZipFile(tmp,'w',ZIP_DEFLATED) as out:
  for info in z.infolist():out.writestr(info,updates.get(info.filename,z.read(info.filename)))
tmp.replace(P)
w=load_workbook(P);v=load_workbook('revisions/submission-ready/Supplementary-Table-8.xlsx')
for a,b in [('NNRTI tyrosines','Table 8'),('Tyrosine atoms','Atoms')]:assert list(w[a].values)==list(v[b].values)
print('Only new-sheet header wrapping, row height and column widths adjusted; all values unchanged.')
