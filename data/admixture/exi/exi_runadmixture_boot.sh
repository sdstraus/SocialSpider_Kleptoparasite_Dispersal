for K in 4; \
do ./admixture \
-B \
--cv=10 \
exi.scaffolds.bed \
$K | tee log_boot${K}.out; \
done
