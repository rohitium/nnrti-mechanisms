from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from lxml import etree as E
from copy import deepcopy
from openpyxl import load_workbook
S=Path('revisions/submission-ready');O=Path('/tmp/nnrti-nov4-stage');ns='http://schemas.openxmlformats.org/spreadsheetml/2006/main';rel='http://schemas.openxmlformats.org/officeDocument/2006/relationships';pkg='http://schemas.openxmlformats.org/package/2006/relationships'
authored=load_workbook(O/'table7-authored.xlsx',read_only=True,data_only=False);old8=load_workbook(S/'Supplementary-Table-8.xlsx',read_only=True,data_only=False)
for name,old in [('NNRTI tyrosines','Table 8'),('Tyrosine atoms','Atoms')]:assert [row[:old8[old].max_column] for row in authored[name].values]==list(old8[old].values)
with ZipFile(S/'Supplementary-Table-7.xlsx') as a,ZipFile(S/'Supplementary-Table-8.xlsx') as b:
 assert a.read('xl/styles.xml')==b.read('xl/styles.xml')
 book=E.fromstring(a.read('xl/workbook.xml'));sheets=book.find('{'+ns+'}sheets');rels=E.fromstring(a.read('xl/_rels/workbook.xml.rels'));ct=E.fromstring(a.read('[Content_Types].xml'));updates={}
 for i,name in enumerate(['NNRTI tyrosines','Tyrosine atoms'],7):
  rid='rId'+str(2+i);node=E.SubElement(sheets,'{'+ns+'}sheet');node.set('name',name);node.set('sheetId',str(i));node.set('{'+rel+'}id',rid)
  node=E.SubElement(rels,'{'+pkg+'}Relationship');node.set('Type',rel+'/worksheet');node.set('Target','/xl/worksheets/sheet'+str(i)+'.xml');node.set('Id',rid)
  node=E.SubElement(ct,'{http://schemas.openxmlformats.org/package/2006/content-types}Override');node.set('PartName','/xl/worksheets/sheet'+str(i)+'.xml');node.set('ContentType','application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml')
  updates['xl/worksheets/sheet'+str(i)+'.xml']=b.read('xl/worksheets/sheet'+str(i-6)+'.xml')
 for key,root in [('xl/workbook.xml',book),('xl/_rels/workbook.xml.rels',rels),('[Content_Types].xml',ct)]:updates[key]=E.tostring(root)
 with ZipFile(O/'Supplementary-Table-7.xlsx','w',ZIP_DEFLATED) as out:
  for info in a.infolist():out.writestr(info,updates.pop(info.filename,a.read(info.filename)))
  for name,data in updates.items():out.writestr(name,data)
w=load_workbook(O/'Supplementary-Table-7.xlsx');old7=load_workbook(S/'Supplementary-Table-7.xlsx');old8=load_workbook(S/'Supplementary-Table-8.xlsx')
for x,y in [(sh,old7[sh.title]) for sh in list(w)[:6]]+[(w['NNRTI tyrosines'],old8['Table 8']),(w['Tyrosine atoms'],old8['Atoms'])]:
 assert list(x.values)==list(y.values)
 assert x.freeze_panes==y.freeze_panes and x.auto_filter.ref==y.auto_filter.ref
 assert {k:v.width for k,v in x.column_dimensions.items()}=={k:v.width for k,v in y.column_dimensions.items()}
 for row in y:
  for cell in row:assert x[cell.coordinate]._style==cell._style,(x.title,cell.coordinate)
print('Merged 8 sheets; every source value, cell style, column width, freeze pane and filter preserved.')
