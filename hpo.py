import re,collections
def load_obo(p):
    terms={};cur=None
    for line in open(p):
        line=line.strip()
        if line=="[Term]": cur={"is_a":[], "obs":False}; continue
        if line.startswith("[") : cur=None; continue
        if cur is None: continue
        if line.startswith("id: "): cur["id"]=line[4:]; terms[cur["id"]]=cur
        elif line.startswith("name: "): cur["name"]=line[6:]
        elif line.startswith("is_a: "): cur["is_a"].append(line[6:16])
        elif line.startswith("is_obsolete: true"): cur["obs"]=True
        elif line.startswith("alt_id: "): cur.setdefault("alt",[]).append(line[8:])
    children=collections.defaultdict(set)
    for t,d in terms.items():
        for p_ in d["is_a"]: children[p_].add(t)
    alt={a:t for t,d in terms.items() for a in d.get("alt",[])}
    return terms,children,alt
def desc(root,children):
    out={root};st=[root]
    while st:
        x=st.pop()
        for c in children.get(x,()):
            if c not in out: out.add(c);st.append(c)
    return out
