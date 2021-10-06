ligo-skymap-stats \
`# Write output to bayestar.tsv.` \
-o bayestar.tsv \
`# Include this option to enable P-P plots.` \
--database coinc.sqlite \
`# Read all sky maps in this directory.` \
*.fits \
`# Optional: calculate the 50% and 90% credible areas.` \
--contour 50 90 \
`# Optional: calculate the probability contained within the smallest` \
`# credible regions of 10 and 100 deg2.` \
--area 10 100 \
`# Optional: count the number of disjoint patches on the sky.` \
`# WARNING: this option makes the script very slow!` \
--modes \
`# Optional, but highly recommended: analyze sky maps using multiple` \
`# threads. In this example, we use 8 worker processes.` \
-j 4
