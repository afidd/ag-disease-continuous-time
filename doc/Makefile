ADCT=../src/disease
NAADSMTRACE=/work/ajd27/naadsmdata/HPAI_Kershaw_Parameters_no_controls_airborne_only_p5_runs1000_2.out
NAADSMDATA=naadsm.h5
ADCTDATA=run.h5
SCENARIO=../data/HPAI_Kershaw_Parameters_no_controls_airborne_only_p5.xml
HERDFILE=../data/HPAI_Kershaw_Unit_no_controls.xml

.PHONY: ALL rpackages

ALL: clinical.pdf latent.pdf susceptible.pdf total_outbreak_compare.pdf

clinical.pdf: clinical.csv
	R --no-save --args "clinical" < individual_dist.R

latent.pdf: latent.csv
	R --no-save --args "latent" < individual_dist.R

susceptible.pdf: susceptible.csv
	R --no-save --args "susceptible" < individual_dist.R

total_outbreak_compare.pdf: adct_hist.csv naadsm_hist.csv
	R --no-save --args adct_hist.csv naadsm_hist.csv < survive_totals.R

clinical.csv latent.csv susceptible.csv: $(NAADSMDATA)
	python residence_discrete.py --input $(NAADSMDATA)

naadsm_hist.csv: $(NAADSMDATA)
	python outbreaksize.py --input $(NAADSMDATA) --output naadsm_hist.csv

adct_hist.csv: $(ADCTDATA)
	python outbreaksize.py --input $(ADCTDATA) --output adct_hist.csv

$(NAADSMDATA): $(NAADSMTRACE)
	python read_naadsm.py --multiple --input $(NAADSMTRACE) --output $(NAADSMDATA)

$(ADCTDATA): $(ADCT)
	$(ADCT) -s $(SCENARIO)  -h $(HERDFILE) --seed 332011 -j 4 --runcnt 1000 --loglevel=warning

# This target helps to initialize a new R installation
# with all of the libraries used in the R scripts
rpackages:
	for p in `grep -h library *.R | sort | uniq | cut -c9- | cut -d')' -f1`; \
	do                                                                       \
		R -e "install.packages(\"$$p\", repos=\"http://lib.stat.cmu.edu/R/CRAN/\")"; \
	done;
