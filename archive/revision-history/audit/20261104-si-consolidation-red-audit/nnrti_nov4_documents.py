from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from lxml import etree as E
from copy import deepcopy
from difflib import SequenceMatcher
from collections import Counter
from hashlib import sha256
import re,json,bisect
ROOT=Path('/Users/rohitpro/Career/00_Github/nnrti-mechanisms');S=ROOT/'revisions/submission-ready';O=Path('/tmp/nnrti-nov4-stage');O.mkdir(exist_ok=True)
W='http://schemas.openxmlformats.org/wordprocessingml/2006/main';M='http://schemas.openxmlformats.org/officeDocument/2006/math';NS={'w':W,'m':M};T=lambda x:'{'+W+'}'+x
XMLSPACE='{http://www.w3.org/XML/1998/namespace}space'
def load(path):
 with ZipFile(path) as z:return E.fromstring(z.read('word/document.xml'))
def txt(p):return ''.join(p.xpath('.//w:t/text()',namespaces=NS))
def norm(s):return re.sub(r'\s+',' ',s).strip()
def paras(root):return root.findall('.//w:body//w:p',NS)
def find(root,s):
 hits=[p for p in paras(root) if s in txt(p)];assert len(hits)==1,(s,len(hits));return hits[0]
def color(r,val):
 pr=r.find(T('rPr'))
 if pr is None:pr=E.Element(T('rPr'));r.insert(0,pr)
 c=pr.find(T('color'))
 if c is None:c=E.SubElement(pr,T('color'))
 c.attrib.clear();c.set(T('val'),val)
def replace(p,old,new):
 nodes=p.xpath('.//w:t',namespaces=NS);full=''.join(n.text or '' for n in nodes);assert old in full,(old,full[:70])
 # All target spans are plain text, outside citation instruction fields.
 pos=full.index(old);end=pos+len(old);offset=0;hit=[]
 for n in nodes:
  a,b=offset,offset+len(n.text or '');offset=b
  if a<end and b>pos:hit.append((n,a,b))
 assert hit
 for i,(n,a,b) in enumerate(hit):
  x=n.text or '';pre=x[:max(0,pos-a)];post=x[min(len(x),end-a):]
  n.text=pre+(new if i==0 else '')+post;n.set(XMLSPACE,'preserve')
def rewrite(p,s):
 assert not p.xpath('.//w:instrText|.//w:drawing',namespaces=NS)
 template=next(iter(p.findall(T('r'))),None);pr=deepcopy(template.find(T('rPr'))) if template is not None and template.find(T('rPr')) is not None else None
 for ch in list(p):
  if ch.tag!=T('pPr'):p.remove(ch)
 r=E.SubElement(p,T('r'))
 if pr is not None:r.append(pr)
 n=E.SubElement(r,T('t'));n.text=s;n.set(XMLSPACE,'preserve');color(r,'000000')
def fieldgroups(p):
 groups=[];current=[];depth=0
 for child in p:
  flags=child.xpath('.//w:fldChar/@w:fldCharType',namespaces=NS)
  if 'begin' in flags:depth+=1
  if depth:current.append(child)
  if 'end' in flags:
   depth-=1
   if depth==0:
    code=''.join(x for ch in current for x in ch.xpath('.//w:instrText/text()',namespaces=NS))
    if 'CSL_CITATION' in code:
     try:obj=json.loads(code[code.index('{'):]);groups.append((obj,list(current)))
     except Exception as e:raise ValueError(code[:120]) from e
    current=[]
 return groups
def itemkey(item):
 d=item.get('itemData',{});return norm((d.get('DOI') or d.get('title') or str(item.get('uris') or item.get('id')))).lower()
ms=load(S/'Manuscript-Marked.docx');rp=load(S/'Response-to-Reviewers.docx');si=load(S/'Supporting-Information-Clean.docx');cover=load(S/'Supporting-Information-Cover-Page-Clean.docx')
original=load(ROOT/'paper/submission/DorDRM-MD-09-02-26-ACS.docx')
initial={n:load(S/n) for n in ['Manuscript-Marked.docx','Response-to-Reviewers.docx','Supporting-Information-Clean.docx','Supporting-Information-Cover-Page-Clean.docx']}
# Preserve the author's exact new SI wording from the response letter.
for start,oldstart in [('We address this question in two steps:','We address this question in two steps:'),('To place the agreement between the two DOR-bound structures','To calibrate these alignment results')]:
 s=txt(find(rp,start)).strip().strip('“”');rewrite(find(si,oldstart),s)
for root in [ms,rp]:
 p=find(root,'relax into the free energy minimum corresponding to the WT RT-DOR simulations performed in this study')
 replace(p,'relax into the free energy minimum corresponding to the WT RT-DOR simulations performed in this study','return to the sampled WT RT–DOR binding-pose basin')
# All citations/inventories now refer to one multi-sheet Table 7.
for root in [ms,rp,si,cover]:
 for p in list(paras(root)):
  for old,new in [('Supplementary Tables 7 and 8','Supplementary Table 7'),('Supplementary Tables 1–8','Supplementary Tables 1–7'),('Tables 7–8 contain the independent-structure and tyrosine comparisons','Table 7 contains the independent-structure and tyrosine comparisons on separate sheets')]:
   if old in txt(p):replace(p,old,new)
# Consolidate the publication inventories, keeping only Table 7.
p=find(ms,'Supplementary Table 7: Comparison of');rewrite(p,'Supplementary Table 7: Comparison of the doravirine-bound NNIBP between PDB 4NCG and 7Z2G, and Tyr181/Tyr188 comparisons across RT–NNRTI structures, on separate sheets (XLSX)')
p=find(ms,'Supplementary Table 8:');p.getparent().remove(p)
p=find(cover,'Supplementary Table 7 (XLSX).');rewrite(p,'Supplementary Table 7 (XLSX). Comparison of the doravirine-bound NNIBP between PDB 4NCG and 7Z2G, with separate sheets comparing Tyr181 and Tyr188 across RT–NNRTI structures.')
p=find(cover,'Supplementary Table 8 (XLSX).');p.getparent().remove(p)
# Synchronize existing native response citation fields with the author's refreshed manuscript.
mainfields={obj['citationID']:group for p in paras(ms) for obj,group in fieldgroups(p)}
for p in paras(rp):
 for obj,group in fieldgroups(p):
  if obj['citationID'] not in mainfields:continue
  pos=p.index(group[0]);italic=any(ch.xpath('.//w:i[@w:val="1"]',namespaces=NS) for ch in group)
  for ch in group:p.remove(ch)
  for j,ch in enumerate(mainfields[obj['citationID']]):
   rr=deepcopy(ch)
   for r in rr.iter(T('r')):
    color(r,'000000');pr=r.find(T('rPr'))
    for tag,value in [('b','0'),('i','1' if italic else '0')]:
     v=pr.find(T(tag))
     if v is None:v=E.SubElement(pr,T(tag))
     v.set(T('val'),value)
   p.insert(pos+j,rr)
# Black all response prose, retaining the author's bold/italic distinctions.
for r in rp.xpath('//w:r|//m:r',namespaces=NS):color(r,'000000')
# Citation-aware, character-precise red font comparison.
def citation_runs(root):
 excluded=set();groups=[]
 for p in paras(root):
  for obj,group in fieldgroups(p):
   runs=[r for ch in group for r in ([ch] if ch.tag==T('r') else ch.iter(T('r')))]
   excluded.update(runs);groups.append((obj,runs))
 return excluded,groups
oldexcluded,oldgroups=citation_runs(original);newexcluded,newgroups=citation_runs(ms)
olditems={itemkey(item) for obj,_ in oldgroups for item in obj.get('citationItems',[])}
def stream(root,excluded):
 chunks=[];locations=[];tokens=[];offset=0;paragraphs=[]
 for p in paras(root):
  begin=offset
  for n in p.xpath('.//w:t|.//m:t',namespaces=NS):
   if n.getparent() in excluded:continue
   s=n.text or '';chunks.append(s);locations.append((offset,offset+len(s),n));offset+=len(s)
  chunks.append('\n');offset+=1;paragraphs.append((begin,offset,p))
 full=''.join(chunks);matches=list(re.finditer(r'\w+|[^\w\s]',full));values=[m.group() for m in matches]
 # Figure/table/section renumbering is a cross-reference update, not added prose.
 for m in re.finditer(r'\b(?:Supplementary\s+)?(?:Figures?|Tables?|Sections?)\s+\d+[A-Za-z]?(?:\s*(?:[–,]|and)\s*\d+[A-Za-z]?)*',full):
  for k,t in enumerate(matches):
   if t.start()>=m.end():break
   if t.start()>=m.start() and re.fullmatch(r'\d+[A-Za-z]?',t.group()):values[k]='#REF'
 return full,values,matches,locations,paragraphs
oldfull,a,am,al,ap=stream(original,oldexcluded);newfull,b,bm,bl,bp=stream(ms,newexcluded)
added=[]
for op,i,j,k,l in SequenceMatcher(None,a,b,autojunk=False).get_opcodes():
 if op in ['replace','insert']:
  added.extend((m.start(),m.end()) for m in bm[k:l])
mask=bytearray(len(newfull))
for x,y in added:mask[x:y]=b'\1'*(y-x)
# Exempt unchanged/moved complete paragraphs, ignoring bibliography labels and reference caches.
def pkey(s):
 s=re.sub(r'^\(\d+\)\s*','',norm(s));s=re.sub(r'\b((?:Supplementary\s+)?(?:Figures?|Tables?|Sections?))\s+\d+[A-Za-z]?(?:\s*(?:[–,]|and)\s*\d+[A-Za-z]?)*',r'\1 #REF',s)
 return s
oldkeys={pkey(oldfull[x:y]) for x,y,p in ap}
for x,y,p in bp:
 if pkey(newfull[x:y]) in oldkeys:mask[x:y]=b'\0'*(y-x)
# Compare each changed paragraph with its closest original paragraph. A global
# sequence alone can misalign when tables or paragraphs move around new figures.
oldstarts=[m.start() for m in am];newstarts=[m.start() for m in bm]
oldsegments=[]
for x,y,p in ap:
 i,j=bisect.bisect_left(oldstarts,x),bisect.bisect_left(oldstarts,y)
 if i<j:oldsegments.append(a[i:j])
oldngrams={tuple(at[k:k+8]) for at in oldsegments for k in range(max(0,len(at)-7))}
for x,y,p in bp:
 if pkey(newfull[x:y]) in oldkeys:continue
 i,j=bisect.bisect_left(newstarts,x),bisect.bisect_left(newstarts,y)
 bt=b[i:j]
 if not bt:continue
 candidates=sorted(((SequenceMatcher(None,at,bt,autojunk=False).quick_ratio(),at) for at in oldsegments),key=lambda v:v[0],reverse=True)[:8]
 ranked=[(SequenceMatcher(None,at,bt,autojunk=False).ratio(),at) for _,at in candidates]
 score,at=max(ranked,key=lambda v:v[0])
 mask[x:y]=b'\1'*(y-x)
 if score>=.40:
  for block in SequenceMatcher(None,at,bt,autojunk=False).get_matching_blocks():
   for k in range(i+block.b,i+block.b+block.size):
    m=bm[k];mask[m.start():m.end()]=b'\0'*(m.end()-m.start())
 # Recognize unchanged phrases even when clauses change order within a paragraph.
 for z in range(max(0,len(bt)-7)):
  if tuple(bt[z:z+8]) in oldngrams:
   for k in range(i+z,i+z+8):
    m=bm[k];mask[m.start():m.end()]=b'\0'*(m.end()-m.start())

# Red whitespace only when it joins two red tokens; this avoids isolated color fragments.
for x,y,n in bl:
 s=newfull[x:y]
 for m in re.finditer(r'\s+',s):
  u,v=x+m.start(),x+m.end()
  mask[u:v]=b'\0'*(v-u)
  if u>0 and v<len(mask) and mask[u-1] and mask[v]:mask[u:v]=b'\1'*(v-u)
for r in ms.xpath('//w:r|//m:r',namespaces=NS):color(r,'000000')
redreport=[];fallback=[]
for x,y,n in bl:
 if x==y:continue
 r=n.getparent();bits=mask[x:y];s=n.text or ''
 # Split only ordinary text runs; all embedded objects and live fields remain intact.
 allowed={T('rPr'),T('t'),T('lastRenderedPageBreak'),'{'+M+'}rPr','{'+M+'}t'}
 if len(set(bits))>1 and all(ch.tag in allowed for ch in r) and len(r.xpath('./w:t|./m:t',namespaces=NS))==1:
  parent=r.getparent();idx=parent.index(r);parts=[];start=0
  for k in range(1,len(bits)+1):
   if k==len(bits) or bits[k]!=bits[start]:
    rr=deepcopy(r);nn=rr.find(n.tag);nn.text=s[start:k];nn.set(XMLSPACE,'preserve');color(rr,'FF0000' if bits[start] else '000000');parts.append(rr);start=k
  parent.remove(r)
  for k,rr in enumerate(parts):
   if k:
    for hint in rr.findall(T('lastRenderedPageBreak')):rr.remove(hint)
   parent.insert(idx+k,rr)
 else:
  if len(set(bits))>1 and not re.match(r'^\(\d+\)',txt(r.getparent()).strip()):fallback.append(s[:80])
  color(r,'FF0000' if all(bits) else '000000')
for obj,runs in newgroups:
 isnew=any(itemkey(item) not in olditems for item in obj.get('citationItems',[]))
 for r in runs:color(r,'FF0000' if isnew else '000000')
# Existing bibliography entries may move as a block; all of their text stays black.
oldbib={pkey(txt(p)) for p in paras(original) if re.match(r'^\(\d+\)',txt(p).strip())}
for p in paras(ms):
 if re.match(r'^\(\d+\)',txt(p).strip()):
  for r in p.xpath('.//w:r',namespaces=NS):color(r,'000000' if pkey(txt(p)) in oldbib else 'FF0000')
for p in paras(ms):
 red=''.join(n.text or '' for n in p.xpath('.//w:r[w:rPr/w:color[@w:val="FF0000"]]/w:t',namespaces=NS))
 if red:redreport.append({'paragraph':txt(p)[:90],'red_text':red})
assert not fallback,fallback
# Regression guard: the author's example is wholly original prose apart from renumbered references.
p=find(ms,'For each system, we ran three independent 100ns')
start=txt(p).index('For each system,');nodes=p.xpath('.//w:t',namespaces=NS);offset=0
for n in nodes:
 end=offset+len(n.text or '')
 if end>start:assert n.getparent().find(T('rPr')).find(T('color')).get(T('val'))=='000000',(n.text,'should be unchanged')
 offset=end
# Original paragraph text, equations and all Zotero instructions are untouched by recoloring.
assert Counter(initial['Manuscript-Marked.docx'].xpath('//w:instrText/text()',namespaces=NS))==Counter(ms.xpath('//w:instrText/text()',namespaces=NS))
# Response quotes must agree with the revised sources, including the author-edited new SI text.
joined=norm(' '.join(txt(p) for root in [ms,si] for p in paras(root)))
qfail=[]
for p in paras(rp):
 s=txt(p).strip()
 if s.startswith('“') or (s.endswith('”') and '“' not in s):
  if s.startswith(('Manuscript,','Supplementary Text,')):continue
  v=s.strip('“”').replace('…','').strip()
  if v and norm(v) not in joined:qfail.append(v[:180])
assert not qfail,qfail
for name,root in [('Manuscript-Marked.docx',ms),('Response-to-Reviewers.docx',rp),('Supporting-Information-Clean.docx',si),('Supporting-Information-Cover-Page-Clean.docx',cover)]:
 assert not re.search(r'Supplementary Tables? 8\b|Supplementary Tables 7 and 8|Tables 7–8', ' '.join(txt(p) for p in paras(root))),name
 def save(dest,clean=False):
  updates={'word/document.xml':E.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)}
  with ZipFile(S/name) as zi,ZipFile(dest,'w',ZIP_DEFLATED) as zo:
   for info in zi.infolist():
    data=updates.get(info.filename,zi.read(info.filename))
    if clean and info.filename.startswith('word/') and info.filename.endswith('.xml'):
     try:r=E.fromstring(data)
     except E.XMLSyntaxError:zo.writestr(info,data);continue
     for c in r.findall('.//w:color',NS):c.attrib.clear();c.set(T('val'),'000000')
     for h in r.xpath('//w:highlight|//w:trackRevisions|//w:rPrChange|//w:pPrChange',namespaces=NS):h.getparent().remove(h)
     assert not r.xpath('//w:ins|//w:del|//w:moveTo|//w:moveFrom|//w:commentReference',namespaces=NS)
     data=E.tostring(r,xml_declaration=True,encoding='UTF-8',standalone=True)
    zo.writestr(info,data)
 save(O/name,clean='Clean' in name)
 if name=='Manuscript-Marked.docx':save(O/'Manuscript-Clean.docx',clean=True)
(O/'red-font-audit.json').write_text(json.dumps({'red_paragraphs':redreport,'old_citation_items':len(olditems),'citation_groups':len(newgroups),'mixed_run_fallbacks':fallback,'unchanged_example_verified_black':True},indent=2,ensure_ascii=False))
print('Staged five DOCX files; all response quotations agree; example paragraph is black; reference fields preserved.')
print('Paragraphs with red text:',len(redreport))
