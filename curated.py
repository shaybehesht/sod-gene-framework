"""Curated HPO terms for the three SOD features (Table S3).

HP:0100842 (Septo-optic dysplasia) is listed under midline and optic for the earlier version of
the framework. With SOD_SEP=1 (the published analysis) build_sets.py removes it from both, and it
is used only as a separate inclusion flag. Published counts: midline 12, ONH 4, pituitary 19.
"""
MID=['HP:0001273', 'HP:0001274', 'HP:0001331', 'HP:0001338', 'HP:0002079', 'HP:0006849', 'HP:0006989', 'HP:0007370', 'HP:0007375', 'HP:0030283', 'HP:0033725', 'HP:0100842', 'HP:0200012']
OPT=['HP:0000609', 'HP:0007766', 'HP:0008058', 'HP:0012521', 'HP:0100842']
PIT=['HP:0000871', 'HP:0010625', 'HP:0010626', 'HP:0010627', 'HP:0011747', 'HP:0011751', 'HP:0011753', 'HP:0011755', 'HP:0011756', 'HP:0011757', 'HP:0012503', 'HP:0012504', 'HP:0012506', 'HP:0012731', 'HP:0034976', 'HP:0034977', 'HP:0034978', 'HP:0040075', 'HP:0410279']
MID=set(MID);OPT=set(OPT);PIT=set(PIT)
