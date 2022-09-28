for K in 1 2 3 4 5 6 7 8 9 10 11 12; \
do ./admixture \
--cv=10 \
exi_subset.js.bed \
$K | tee log${K}.out; \
done
