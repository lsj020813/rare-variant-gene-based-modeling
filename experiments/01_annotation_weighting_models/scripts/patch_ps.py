p='partial_summary.sh'; s=open(p).read()
old="""    g=d.get('gate_uniform'); 
    if g and g['gate_max_abs_dpip'] not in ('NA','') and float(g['gate_max_abs_dpip'])<1e-6: gate_ok+=1"""
new="""    g=d.get('gate_uniform')
    if g is None:
        gf=f.replace('.arms.tsv','.gate.arms.tsv')
        if os.path.exists(gf): g={r['arm']:r for r in csv.DictReader(open(gf),delimiter='\\t')}.get('gate_uniform')
    if g and g['gate_max_abs_dpip'] not in ('NA','') and float(g['gate_max_abs_dpip'])<1e-6: gate_ok+=1"""
assert old in s; s=s.replace(old,new).replace('import glob, csv','import glob, csv, os'); open(p,'w').write(s); print('patched')
